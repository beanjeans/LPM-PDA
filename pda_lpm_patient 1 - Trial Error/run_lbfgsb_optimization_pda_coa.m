%% RUN_LBFGSB_OPTIMIZATION_PDA_COA
% =========================================================================
% L-BFGS-B STYLE BOUNDED OPTIMIZATION — PDA-CoA LPM
%
% PURPOSE:
%   Calibrate influential model parameters (identified from Sobol GSA)
%   by minimizing a weighted normalized least-squares error between
%   simulated and clinical haemodynamic measurements.
%
%   Expanded objective: MAP, SV, SBP, DBP, PDA gradient and, when
%   measured, CoA gradient, Qp/Qs, LV EF, and LV EDV.
%
%   Uses MATLAB fmincon with the interior-point algorithm, bound constraints,
%   normalized variables, and a limited-memory BFGS Hessian approximation.
%
%   Patient 1 is PDA-only. MAP, SV, SBP, DBP, and dP_PDA are calibrated.
%   Virtual CoA geometry remains fixed near zero during calibration and is
%   introduced later during the virtual CoA experiment.
%
% PARAMETERS OPTIMIZED:
%   Emax_lv, R_systemic, R_shunt_pda, R_pa, C_sys, and C_ao.
% WORKFLOW:
%   1. Load patient clinical data (non-interactive, batch mode)
%   2. Build baseline model parameters
%   3. Run baseline simulation & compute baseline error
%   4. Run bounded quasi-Newton optimization (fmincon)
%   5. Apply optimized parameters & run final simulation
%   6. Save results to CSV
%   7. Generate plots
%   8. Report all active clinical calibration targets
%
% OUTPUT FILES (saved to results/optimization/):
%   optimized_parameters.csv    — optimized vs baseline parameter values
%   objective_history.csv       — J and x per optimizer iteration
%   opt_convergence.png         — objective convergence + parameter traces
%   opt_clinical_comparison.png — MAP and SV before vs after optimization
%
% IMPORTANT:
%   This script does NOT modify any existing model files.
%   It calls: build_coa_params, system_rhs_pda_coa,
%             compute_clinical_indices, default_parameters
%
% REFERENCES:
%   [1] Nocedal & Wright (2006). Numerical Optimization. Springer.
%   [2] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%
% AUTHOR:   Optimization Extension — Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% =========================================================================

clear; clc; close all;

%% Add all required paths
addpath('config', 'models', 'solvers', 'utils', 'tests', 'gsa', 'optimization');

fprintf('=================================================================\n');
fprintf('   L-BFGS-B BOUNDED PARAMETER OPTIMIZATION\n');
fprintf('   PDA-CoA Lumped Parameter Cardiovascular Model\n');
fprintf('=================================================================\n\n');

% =========================================================================
%  SECTION A — USER CONFIGURATION  (edit here only)
% =========================================================================

%% A1. Patient selection
csv_path    = fullfile('config', 'patient_data.csv');
patient_idx = 1;    % Patient row index in patient_data.csv (1-based)

%% A2. Parameters to optimize (Expanded set to prevent parameter non-identifiability)
%  Includes baseline top GSA rankers plus parallel shunt and compliance pathways.
opt_param_names = {
    'Emax_lv'         % LV peak elastance         [mmHg/mL]    — Systolic contractility
    'R_systemic'      % Total systemic resistance [mmHg·s/mL]  — Downstream SVR
    'R_shunt_pda'     % PDA shunt resistance      [mmHg·s/mL]  — Trans-ductal flow control
    'R_pa'            % Pulmonary resistance      [mmHg·s/mL]  — Downstream PVR
    'C_sys'           % Systemic compliance       [mL/mmHg]    — Peripheral storage
    'C_ao'            % Aortic compliance         [mL/mmHg]    — Proximal Windkessel
};

%% A3. Parameter bounds  [lower, upper]
%  Physiological neonatal bounds (matching allometric scaling ranges).
opt_bounds = [
%   Lower    Upper
    3.0,     20.0   % Emax_lv      [mmHg/mL]
    1.0,     15.0   % R_systemic   [mmHg·s/mL]
    0.5,     20.0   % R_shunt_pda  [mmHg·s/mL]
    0.05,    1.0    % R_pa         [mmHg·s/mL]
    0.05,    0.5    % C_sys        [mL/mmHg]
    0.0001,  0.002  % C_ao         [mL/mmHg]
];

%% A4. Objective mode and weights
objective_mode = 'expanded_physiological_targets';
primary_optimization_metrics = {'MAP', 'SV', 'SBP', 'DBP', 'dP_PDA', ...
    'dP_CoA_peak', 'Qp_Qs', 'EF_lv', 'EDV_lv'};

weights.MAP          = 10.0;  % Primary systemic pressure target
weights.SV           = 6.0;   % Primary stroke volume target
weights.SBP          = 4.0;   % Constrains systolic peak and elastance
weights.DBP          = 4.0;   % Constrains diastolic decay and compliance
weights.PP           = 0.0;   % Avoid double-counting SBP and DBP
weights.dP_PDA       = 3.0;   % Ductal pressure gradient (prevents runaway shunt)
weights.dP_CoA_peak  = 0.0;   % P01 has no measured CoA gradient
weights.dP_CoA_mean  = 0.0;   % Disabled
weights.Q_CoA_frac   = 0.0;   % Disabled
weights.Qp_Qs        = 0.0;   % Enable only when clinically measured
weights.EF_lv        = 0.0;   % Enable only when clinically measured
weights.EDV_lv       = 0.0;   % Enable only when clinically measured
weights.parameter_prior = 0.1; % Weak regularization for identifiability

%% A5. Fixed CoA geometry (used if stenosis_pct / coa_length_mm NOT in opt_param_names)
default_stenosis_pct   =  0.1;   % [%]  — starting geometry
default_coa_length_mm  =  5.0;   % [mm] — mid-range scenario

%% A6. Solver settings
n_warmup = 15;    % Longer settling period for a smoother objective
n_report = 3;     % Reporting cycles used to calculate outputs
penalty  = 1e6;   % Objective value returned on ODE/build failure

%% A7. fmincon options
% fmincon optimizes normalized variables z in [0,1].
fmincon_opts = optimoptions('fmincon', ...
    'Algorithm',              'interior-point', ...
    'HessianApproximation',   'lbfgs', ...
    'Display',                'iter', ...
    'MaxIterations',          300, ...
    'MaxFunctionEvaluations', 5000, ...
    'OptimalityTolerance',    1e-4, ...
    'StepTolerance',          1e-6, ...
    'FiniteDifferenceType',   'central', ...
    'FiniteDifferenceStepSize', 1e-3, ...
    'ScaleProblem',           true, ...
    'OutputFcn',              @optimization_output_callback);

%% A8. Output directory
results_dir = fullfile('results', 'optimization');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% =========================================================================
%  STEP 1 — LOAD PATIENT DATA (non-interactive)
% =========================================================================
fprintf('STEP 1: Loading patient data...\n');

patient_table = readtable(csv_path);
row = patient_table(patient_idx, :);

clinical.patient_id         = row.PatientID{1};
clinical.age_days           = row.Age(1);
clinical.sex                = row.Sex{1};
clinical.height_cm          = row.TB(1);
clinical.weight_g           = row.BB(1);
clinical.BSA_m2             = row.BSA(1);
clinical.HR_bpm             = row.HeartRate(1);
clinical.SV_mL              = row.StrokeVolume(1);
clinical.P_ao_sys_mmHg      = row.SSAP(1);
clinical.P_ao_dia_mmHg      = row.SDAP(1);
clinical.P_ao_mean_mmHg     = row.MAP(1);
clinical.D_shunt_pda_mm     = row.DPDA(1);
clinical.D_coa_mm           = row.DCoA(1);
clinical.D_aao_mm           = row.DAAo(1);
clinical.D_dta_mm           = row.DDTA(1);
clinical.D_isthmus_mm       = row.DIsthmus(1);
clinical.D_dao_mm           = row.DDAo(1);
clinical.D_aov_mm           = row.DAoV(1);
clinical.D_pv_mm            = row.DPV(1);
clinical.v_aov_ms           = row.vAoV(1);
clinical.v_pv_psax_ms       = row.vPVpsax(1);
clinical.v_pv_supra_ms      = row.vPVsupra(1);
clinical.v_pda_ms           = row.vPDA(1);
clinical.v_coa_ms           = row.vCoA(1);
clinical.pda_direction      = row.ArahAliranPDA(1);
clinical.dP_aov_mmHg        = row.dPAoV(1);
clinical.dP_pv_mmHg         = row.dPPV(1);
clinical.dP_pda_mmHg        = row.dPPDA(1);
clinical.dP_coa_mmHg        = row.dPCoA(1);

% Optional measured targets. If these columns are absent, the target is
% NaN and the objective function excludes it automatically.
clinical.Qp_Qs     = table_value_or_nan(row, {'QpQs','Qp_Qs'});
clinical.EF_lv_pct = table_value_or_nan(row, {'LVEF','EF','EF_lv'});
clinical.EDV_lv_mL = table_value_or_nan(row, {'LVEDV','EDV_lv','EDV'});

CO_mLs = (clinical.SV_mL * clinical.HR_bpm) / 60;
clinical.CO_Lmin = CO_mLs * (60 / 1000);
clinical.CO_mLs  = CO_mLs;

if clinical.pda_direction == 1
    clinical.P_pa_est_mmHg = max(clinical.P_ao_mean_mmHg - clinical.dP_pda_mmHg, 5);
else
    clinical.P_pa_est_mmHg = clinical.P_ao_mean_mmHg;
end

fprintf('  Patient: %s | HR: %d bpm | MAP: %.1f mmHg | SV: %.2f mL\n', ...
    clinical.patient_id, clinical.HR_bpm, clinical.P_ao_mean_mmHg, clinical.SV_mL);
fprintf('  No measured CoA target is used in optimization (PDA-only patient).\n\n');

% =========================================================================
%  STEP 2 — BUILD BASELINE PARAMETERS
% =========================================================================
fprintf('STEP 2: Building baseline patient parameters...\n');

% Pass patient BW for allometric scaling (Pennati & Fumero 2000)
BW_neo_kg_opt = clinical.weight_g / 1000;   % [g] → [kg]
params_default = default_parameters(BW_neo_kg_opt);
uc = unit_conversion();

% Calibrate from clinical data — mirrors build_patient_params v2.0 (silent)
params_base = params_default;
params_base.HR_bpm    = clinical.HR_bpm;
params_base.T_cardiac = 60 / clinical.HR_bpm;
params_base.Ts1       = 0.3  * sqrt(params_base.T_cardiac);
params_base.Ts2       = 0.45 * sqrt(params_base.T_cardiac);
params_base.R_systemic = clinical.P_ao_mean_mmHg / clinical.CO_mLs;

% 2b. Pulmonary vascular resistance — clinical override
% Allometric b=-1.00 gives R_pa ≈ 3.55 for 987g; clinical P_pa implies ~0.21.
% Assume Qp/Qs ≈ 1.5 for significant L→R PDA (Rudolph 2001).
Q_pul_est_mLs         = clinical.CO_mLs * 1.5;
P_pv_est_mmHg         = max(clinical.P_pa_est_mmHg - 5, 3);
params_base.R_pa       = max(0.05, min(2.0, ...
    (clinical.P_pa_est_mmHg - P_pv_est_mmHg) / Q_pul_est_mLs));
params_base.R_pv_veins = params_base.R_pa;

% 2c. Systemic venous compliance — clinical override
% Allometric gives C_sys ≈ 0.005 mL/mmHg — starves LV preload (SV ≈ 1.5 mL).
% Clinical estimate: C_sys ≈ SV / (MAP - P_ra_ref),  P_ra_ref = 4 mmHg.
C_sys_clinical          = clinical.SV_mL / max(clinical.P_ao_mean_mmHg - 4, 1);
params_base.C_sys       = max(0.05, min(0.5, C_sys_clinical));

P_lv_target    = clinical.P_ao_mean_mmHg * 1.30;
params_base.Emax_lv = max(3.0, min(20.0, P_lv_target / clinical.SV_mL));
params_base.Emin_lv = params_base.Emax_lv * 0.05;
params_base.Emax_rv = max(1.5, min(12.0, params_base.Emax_lv * 0.5));
params_base.Emin_rv = params_base.Emax_rv * 0.05;

% Doppler-derived PDA resistance (Doppler-only; HP discarded — see build_patient_params.m)
D_pda_m      = clinical.D_shunt_pda_mm * uc.mm_to_m;
A_pda_m2     = pi * (D_pda_m / 2)^2;
Q_pda_vel    = max(A_pda_m2 * clinical.v_pda_ms * uc.m3s_to_mLs, 0.5);
R_pda_Doppl  = clinical.dP_pda_mmHg / Q_pda_vel;
params_base.R_shunt_pda = max(0.5, min(20, R_pda_Doppl));
params_base.P_pa_target_mmHg = clinical.P_pa_est_mmHg;

params_base.X0(params_base.idx.P_ao)  = clinical.P_ao_mean_mmHg;
params_base.X0(params_base.idx.P_sys) = clinical.P_ao_mean_mmHg;
params_base.X0(params_base.idx.P_pa)  = clinical.P_pa_est_mmHg;
params_base.X0(params_base.idx.P_pv)  = P_pv_est_mmHg;
params_base.X0(params_base.idx.P_la)  = max(P_pv_est_mmHg - 1, 2);

fprintf('  R_systemic (baseline):       %.4f mmHg·s/mL\n', params_base.R_systemic);
fprintf('  R_pa       (clin override):  %.4f mmHg·s/mL  (allometric: %.4f)\n', ...
    params_base.R_pa, params_default.R_pa);
fprintf('  C_sys      (clin override):  %.4f mL/mmHg    (allometric: %.4f)\n', ...
    params_base.C_sys, params_default.C_sys);
fprintf('  C_ao       (baseline):       %.6f mL/mmHg\n',   params_base.C_ao);
fprintf('  Emax_lv    (baseline):       %.4f mmHg/mL\n\n', params_base.Emax_lv);

% Diagnostic: verify R_shunt_pda is physiologically meaningful
Q_pda_check = clinical.dP_pda_mmHg / params_base.R_shunt_pda;
fprintf('  R_shunt_pda (Doppler): %.4f mmHg·s/mL\n', params_base.R_shunt_pda);
fprintf('    D_PDA=%.2f mm | v_PDA=%.2f m/s | dP_PDA=%.1f mmHg\n', ...
    clinical.D_shunt_pda_mm, clinical.v_pda_ms, clinical.dP_pda_mmHg);
fprintf('    Q_PDA check: dP/R = %.2f mL/s  |  Doppler A×v = %.2f mL/s\n\n', ...
    Q_pda_check, A_pda_m2 * clinical.v_pda_ms * uc.m3s_to_mLs);

% =========================================================================
%  STEP 3 — BUILD OPTIMIZATION CONFIG STRUCT
% =========================================================================
fprintf('STEP 3: Configuring optimization...\n');

% Extract bounds into separate lb / ub vectors for fmincon
lb = opt_bounds(:, 1);
ub = opt_bounds(:, 2);

assert(size(opt_bounds,1) == numel(opt_param_names), ...
    'Number of parameter bounds must match opt_param_names.');
assert(all(ub > lb), 'Every upper bound must be greater than its lower bound.');

% Build the config struct passed to the objective function
opt_config.param_names   = opt_param_names;
opt_config.fixed_params  = params_base;
opt_config.clinical      = clinical;
opt_config.stenosis_pct  = default_stenosis_pct;
opt_config.coa_length_mm = default_coa_length_mm;
opt_config.weights       = weights;
opt_config.n_warmup      = n_warmup;
opt_config.n_report      = n_report;
opt_config.penalty       = penalty;

% Build initial guess x0 from the baseline params (centre of bounds where unknown)
x0 = zeros(length(opt_param_names), 1);
% Update the switch statement inside Step 3:
for k = 1:length(opt_param_names)
    pname = opt_param_names{k};
    switch pname
        case 'R_systemic',    x0(k) = params_base.R_systemic;
        case 'C_sys',         x0(k) = params_base.C_sys;
        case 'Emax_lv',       x0(k) = params_base.Emax_lv;
        case 'R_shunt_pda',   x0(k) = params_base.R_shunt_pda;
        case 'C_ao',          x0(k) = params_base.C_ao;
        case 'R_pa',          x0(k) = params_base.R_pa;
        otherwise,            x0(k) = (lb(k) + ub(k)) / 2;
    end
    % Clamp x0 to bounds
    x0(k) = max(lb(k), min(ub(k), x0(k)));
end

% Baseline-centred weak prior used to stabilize the expanded parameter set.
opt_config.x_reference = x0;
opt_config.lb = lb;
opt_config.ub = ub;

% Normalize physical parameters x to optimizer variables z in [0,1].
z0 = (x0 - lb) ./ (ub - lb);
to_physical = @(z) lb + z .* (ub - lb);

fprintf('  Optimizing %d parameters:\n', length(opt_param_names));
fprintf('  %-18s  %10s  %8s  %8s\n', 'Parameter', 'x0', 'LB', 'UB');
fprintf('  %s\n', repmat('-', 1, 50));
for k = 1:length(opt_param_names)
    fprintf('  %-18s  %10.4f  %8.4f  %8.4f\n', opt_param_names{k}, x0(k), lb(k), ub(k));
end
fprintf('\n');

% =========================================================================
%  STEP 4 — BASELINE SIMULATION (pre-optimization)
% =========================================================================
fprintf('STEP 4: Running baseline simulation (pre-optimization)...\n');

[J_baseline, baseline_outputs] = objective_lbfgsb_pda_coa(x0, opt_config);
fprintf('  Baseline objective J = %.6f\n\n', J_baseline);

if isempty(baseline_outputs)
    error('RUN_LBFGSB: Baseline simulation failed. Check parameters and patient data.');
end

% =========================================================================
%  STEP 5 — RUN OPTIMIZATION
% =========================================================================
fprintf('STEP 5: Running bounded L-BFGS-B style optimization (fmincon interior-point)...\n');
fprintf('  This will print iteration details below.\n');
fprintf('  fmincon algorithm: interior-point (limited-memory BFGS Hessian)\n\n');

% Shared storage for the output callback
global OPT_HISTORY_X OPT_HISTORY_J OPT_ITER_COUNT;
OPT_HISTORY_X   = z0';         % Store normalized parameter history
OPT_HISTORY_J   = J_baseline;  % 1 × 1
OPT_ITER_COUNT  = 0;

obj_func = @(z) objective_lbfgsb_pda_coa( ...
    to_physical(z), opt_config);

t_opt_start = tic;

[z_opt, J_opt, exitflag, output_struct] = fmincon( ...
    obj_func, z0, ...
    [], [], [], [], ...
    zeros(size(z0)), ones(size(z0)), ...
    [], fmincon_opts);

t_opt_elapsed = toc(t_opt_start);

x_opt = to_physical(z_opt);

fprintf('\n  Optimization finished in %.1f s (%d iterations, %d func evals)\n', ...
    t_opt_elapsed, output_struct.iterations, output_struct.funcCount);
fprintf('  Exit flag: %d  (%s)\n\n', exitflag, exit_flag_message(exitflag));

if isfield(output_struct, 'firstorderopt')
    fprintf('  First-order optimality: %.3e\n', output_struct.firstorderopt);
    if output_struct.firstorderopt > 1e-2
        warning('RUN_LBFGSB:PoorConvergence', ...
            ['Optimization stopped without adequate convergence. ', ...
             'First-order optimality = %.3e.'], output_struct.firstorderopt);
    end
end

% Convert normalized history back to physical units for output and plots.
z_history = OPT_HISTORY_X;
x_history = lb' + z_history .* (ub - lb)';
J_history = OPT_HISTORY_J;

% =========================================================================
%  STEP 6 — POST-OPTIMIZATION SIMULATION
% =========================================================================
fprintf('STEP 6: Running final simulation with optimized parameters...\n');

[J_final, opt_outputs] = objective_lbfgsb_pda_coa(x_opt, opt_config);
fprintf('  Final objective J = %.6f  (improvement: %.2f%%)\n\n', ...
    J_final, 100 * (J_baseline - J_final) / max(J_baseline, eps));

if isempty(opt_outputs)
    warning('RUN_LBFGSB: Post-optimization simulation failed. Reporting x_opt only.');
end

% Also build the optimized params struct for full simulation/reporting
[params_opt, stenosis_opt, coa_length_opt] = apply_optimized_params(x_opt, opt_config);

% =========================================================================
%  STEP 7 — SAVE RESULTS TO CSV
% =========================================================================
fprintf('STEP 7: Saving results...\n');

% --- 7a. Optimized parameters table ---
rows = {};
for k = 1:length(opt_param_names)
    pname = opt_param_names{k};
    rows{end+1, 1} = pname;        %#ok<SAGROW>
    rows{end, 2}   = x0(k);        % Baseline (initial x0)
    rows{end, 3}   = x_opt(k);     % Optimized
    rows{end, 4}   = lb(k);        % Lower bound
    rows{end, 5}   = ub(k);        % Upper bound
    rows{end, 6}   = 100 * (x_opt(k) - x0(k)) / x0(k);  % % change
end

T_params = cell2table(rows, ...
    'VariableNames', {'Parameter', 'Baseline', 'Optimized', 'LB', 'UB', 'ChangePercent'});
params_csv = fullfile(results_dir, 'optimized_parameters.csv');
writetable(T_params, params_csv);
fprintf('  Saved: %s\n', params_csv);

% --- 7b. Objective history ---
n_iter_rec = length(J_history);
iter_nums  = (1:n_iter_rec)';
hist_table_data = [num2cell(iter_nums), num2cell(J_history(:)), num2cell(x_history)];
hist_headers = [{'Iteration', 'J'}, opt_param_names(:)'];
T_history = cell2table(hist_table_data, 'VariableNames', hist_headers);
history_csv = fullfile(results_dir, 'objective_history.csv');
writetable(T_history, history_csv);
fprintf('  Saved: %s\n', history_csv);

% --- 7c. Full workspace ---
% Extract objective breakdowns for convenience
if isstruct(baseline_outputs) && isfield(baseline_outputs, 'objective_breakdown')
    objective_breakdown_baseline = baseline_outputs.objective_breakdown;
else
    objective_breakdown_baseline = struct();
end
if isstruct(opt_outputs) && isfield(opt_outputs, 'objective_breakdown')
    objective_breakdown_final = opt_outputs.objective_breakdown;
else
    objective_breakdown_final = struct();
end

mat_path = fullfile(results_dir, 'optimization_workspace.mat');
save(mat_path, 'x_opt', 'x0', 'z_opt', 'z0', ...
    'J_opt', 'J_baseline', 'J_final', ...
    'x_history', 'z_history', 'J_history', 'opt_config', 'clinical', 'params_opt', ...
    'stenosis_opt', 'coa_length_opt', 'baseline_outputs', 'opt_outputs', ...
    'objective_mode', 'weights', 'primary_optimization_metrics', ...
    'objective_breakdown_baseline', 'objective_breakdown_final');
fprintf('  Saved: %s\n\n', mat_path);

% =========================================================================
%  STEP 8 — GENERATE PLOTS
% =========================================================================
fprintf('STEP 8: Generating plots...\n');

plot_optimization_results(baseline_outputs, opt_outputs, clinical, ...
    x_history, J_history, opt_param_names, results_dir);
end

fprintf('\n');

% =========================================================================
%  STEP 9 — FINAL OPTIMIZATION REPORT
% =========================================================================
fprintf('=========================================================================\n');
fprintf('   PRIMARY OPTIMIZATION OUTPUTS\n');
fprintf('   Patient: %s\n', clinical.patient_id);
fprintf('=========================================================================\n\n');

fprintf('  Objective mode:  %s\n', objective_mode);
fprintf('  Active params:   %s\n', strjoin(opt_param_names, ', '));
fprintf('  Core targets:    MAP, SV, SBP, DBP, dP_PDA\n');
fprintf('  Optional targets are used only when measured and available.\n');
fprintf('\n');

% Helper for safe display
def_nan = @(s, f) get_val_safe(s, f);

% Compute percent errors for direct targets
map_pre  = def_nan(baseline_outputs, 'P_ao_mean');
map_post = def_nan(opt_outputs,      'P_ao_mean');
sv_pre   = def_nan(baseline_outputs, 'SV_lv');
sv_post  = def_nan(opt_outputs,      'SV_lv');
sbp_pre  = def_nan(baseline_outputs, 'P_ao_sys');
sbp_post = def_nan(opt_outputs,      'P_ao_sys');
dbp_pre  = def_nan(baseline_outputs, 'P_ao_dia');
dbp_post = def_nan(opt_outputs,      'P_ao_dia');
pda_pre  = def_nan(baseline_outputs, 'dP_PDA_peak');
pda_post = def_nan(opt_outputs,      'dP_PDA_peak');

map_err_pre  = 100 * (map_pre  - clinical.P_ao_mean_mmHg) / max(abs(clinical.P_ao_mean_mmHg), 1.0);
map_err_post = 100 * (map_post - clinical.P_ao_mean_mmHg) / max(abs(clinical.P_ao_mean_mmHg), 1.0);
sv_err_pre   = 100 * (sv_pre   - clinical.SV_mL)          / max(abs(clinical.SV_mL), 0.1);
sv_err_post  = 100 * (sv_post  - clinical.SV_mL)          / max(abs(clinical.SV_mL), 0.1);
sbp_err_post = 100 * (sbp_post - clinical.P_ao_sys_mmHg)  / max(abs(clinical.P_ao_sys_mmHg), 1.0);
dbp_err_post = 100 * (dbp_post - clinical.P_ao_dia_mmHg)  / max(abs(clinical.P_ao_dia_mmHg), 1.0);
pda_err_post = 100 * (pda_post - clinical.dP_pda_mmHg)    / max(abs(clinical.dP_pda_mmHg), 0.5);

fprintf('  %-10s  %12s  %12s  %12s  %14s\n', ...
    'Target', 'Clinical', 'Pre-Opt', 'Post-Opt', 'Error Post-Opt');
fprintf('  %s\n', repmat('-', 1, 66));
fprintf('  %-10s  %12.2f  %12.2f  %12.2f  %13.2f%%\n', ...
    'MAP (mmHg)', clinical.P_ao_mean_mmHg, map_pre, map_post, map_err_post);
fprintf('  %-10s  %12.2f  %12.2f  %12.2f  %13.2f%%\n', ...
    'SV (mL)',    clinical.SV_mL,          sv_pre,  sv_post,  sv_err_post);
fprintf('  %-10s  %12.2f  %12.2f  %12.2f  %13.2f%%\n', ...
    'SBP', clinical.P_ao_sys_mmHg, sbp_pre, sbp_post, sbp_err_post);
fprintf('  %-10s  %12.2f  %12.2f  %12.2f  %13.2f%%\n', ...
    'DBP', clinical.P_ao_dia_mmHg, dbp_pre, dbp_post, dbp_err_post);
fprintf('  %-10s  %12.2f  %12.2f  %12.2f  %13.2f%%\n', ...
    'dP PDA', clinical.dP_pda_mmHg, pda_pre, pda_post, pda_err_post);
fprintf('  %s\n', repmat('-', 1, 66));
fprintf('\n');
fprintf('  Objective J:    Baseline = %.6f  |  Final = %.6f\n', J_baseline, J_final);
fprintf('  Improvement:    %.2f%%\n', 100*(J_baseline - J_final)/max(J_baseline,eps));
fprintf('\n');
fprintf('  Optimized parameters:\n');
fprintf('    stenosis_pct  = %.2f%%\n', stenosis_opt);
fprintf('    coa_length_mm = %.2f mm  (fixed geometry; not optimized)\n', coa_length_opt);
fprintf('\n');
fprintf('  Expanded physiological targets and weak parameter priors were used.\n');
fprintf('\n');
fprintf('  Results saved to: %s/\n', results_dir);
fprintf('=========================================================================\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function stop = optimization_output_callback(z, optimValues, state)
% Records normalized z and J at each fmincon iteration.
    global OPT_HISTORY_X OPT_HISTORY_J OPT_ITER_COUNT;
    stop = false;
    if strcmp(state, 'iter')
        OPT_ITER_COUNT  = OPT_ITER_COUNT + 1;
        OPT_HISTORY_J   = [OPT_HISTORY_J;   optimValues.fval];
        OPT_HISTORY_X   = [OPT_HISTORY_X;   z(:)'];
    end
end

function msg = exit_flag_message(flag)
% Human-readable fmincon exit flag description
    switch flag
        case  1,  msg = 'Converged to solution';
        case  2,  msg = 'Step size below StepTolerance';
        case  0,  msg = 'Max iterations/evaluations reached';
        case -1,  msg = 'Stopped by output function';
        case -2,  msg = 'No feasible point found';
        otherwise, msg = 'Unknown exit flag';
    end
end

function val = get_val_safe(s, fname)
    if isstruct(s) && isfield(s, fname) && isnumeric(s.(fname))
        val = s.(fname);
    else
        val = NaN;
    end
end

function val = get_str_safe(s, fname)
    if isstruct(s) && isfield(s, fname)
        val = char(s.(fname));
    else
        val = 'N/A';
    end
end

function value = table_value_or_nan(row, candidate_names)
% Return the first finite numeric value from any matching table variable.
value = NaN;
names = row.Properties.VariableNames;
for k = 1:numel(candidate_names)
    idx = find(strcmpi(names, candidate_names{k}), 1);
    if ~isempty(idx)
        raw = row{1, idx};
        if isnumeric(raw) && isscalar(raw) && isfinite(raw)
            value = raw;
            return;
        end
    end
end
end