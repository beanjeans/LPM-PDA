%% RUN_LBFGSB_OPTIMIZATION_PDA_COA
% =========================================================================
% L-BFGS-B STYLE BOUNDED OPTIMIZATION — PDA-CoA LPM
%
% PURPOSE:
%   Calibrate influential model parameters (identified from Sobol GSA)
%   by minimizing a weighted normalized least-squares error between
%   simulated and clinical haemodynamic measurements.
%
%   Objective:
%     J = Σ  w_i × ( (sim_i − clin_i) / clin_i )²
%
%   Uses MATLAB fmincon with the 'interior-point' algorithm and parameter
%   bounds — functionally equivalent to L-BFGS-B bounded optimization.
%   ('interior-point' uses a limited-memory BFGS Hessian internally.)
%
% PARAMETERS OPTIMIZED (from Sobol GSA influential set):
%   C_sys, R_systemic, Emax_lv, stenosis_pct, coa_length_mm
%
% WORKFLOW:
%   1. Load patient clinical data (non-interactive, batch mode)
%   2. Build baseline model parameters
%   3. Run baseline simulation & compute baseline error
%   4. Run bounded quasi-Newton optimization (fmincon)
%   5. Apply optimized parameters & run final simulation
%   6. Save results to CSV
%   7. Generate plots
%   8. Report final CoA clinical outputs
%
% OUTPUT FILES (saved to results/optimization/):
%   optimized_parameters.csv   — optimized vs baseline parameter values
%   objective_history.csv      — J and x per optimizer iteration
%   opt_convergence.png        — objective convergence + parameter traces
%   opt_clinical_comparison.png — before vs after target comparison
%   opt_coa_summary.png        — CoA clinical output summary
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

%% A2. Parameters to optimize (must match struct field names exactly)
%  These are the influential parameters from Sobol GSA.
%  Comment out any parameter you want to fix to its baseline value.
opt_param_names = {
    'R_systemic'      % Total systemic vascular resistance  [mmHg·s/mL]
    'C_sys'           % Systemic venous compliance          [mL/mmHg]
    'Emax_lv'         % LV peak elastance                   [mmHg/mL]
    'stenosis_pct'    % CoA stenosis severity               [%]
};

%% A3. Parameter bounds  [lower, upper]
%  Order must match opt_param_names exactly.
%  Chosen from physiological literature (neonatal ranges).
opt_bounds = [
%   Lower    Upper
    1.0,     20.0    % R_systemic  [mmHg·s/mL]
    0.05,    0.60    % C_sys       [mL/mmHg]
    0.5,     6.0     % Emax_lv     [mmHg/mL]
    5.0,     99.0    % stenosis_pct [%]
];

%% A4. Objective weights
%  Higher weight = this target is more important to match.
%  Set weight to 0 to exclude a target from the objective.
weights.MAP    = 3.0;   % Mean arterial pressure (most reliable clinical target)
weights.SBP    = 2.0;   % Systolic blood pressure
weights.DBP    = 1.5;   % Diastolic blood pressure
weights.SV     = 2.0;   % Stroke volume
weights.dP_PDA = 1.5;   % PDA pressure gradient (Doppler-derived)
weights.dP_CoA = 2.5;   % CoA pressure gradient (Doppler-derived, if available)

%% A5. Fixed CoA geometry (used if stenosis_pct / coa_length_mm NOT in opt_param_names)
default_stenosis_pct   = 50.0;   % [%]  — starting geometry
default_coa_length_mm  =  5.0;   % [mm] — mid-range scenario

%% A6. Solver settings
n_warmup = 6;     % ODE warm-up cycles (validated: steady-state by cycle 6)
n_report = 2;     % ODE reporting cycles
penalty  = 1e6;   % Objective value returned on ODE/build failure

%% A7. fmincon options
% Algorithm: 'interior-point' is the correct bounded L-BFGS-B equivalent in
% MATLAB's fmincon. It uses a limited-memory BFGS Hessian approximation
% internally with explicit bound constraints.
% Valid fmincon algorithms: 'interior-point', 'sqp', 'active-set',
%                           'trust-region-reflective'
fmincon_opts = optimoptions('fmincon', ...
    'Algorithm',              'interior-point', ... % L-BFGS-B equivalent
    'Display',                'iter', ...
    'MaxIterations',          200, ...
    'MaxFunctionEvaluations', 2000, ...
    'OptimalityTolerance',    1e-6, ...
    'StepTolerance',          1e-8, ...
    'FiniteDifferenceType',   'central', ...        % More accurate gradient
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
fprintf('  Clinical CoA gradient (echo): %.1f mmHg\n\n', clinical.dP_coa_mmHg);

% =========================================================================
%  STEP 2 — BUILD BASELINE PARAMETERS
% =========================================================================
fprintf('STEP 2: Building baseline patient parameters...\n');

params_default = default_parameters();
uc = unit_conversion();

% Calibrate from clinical data (mirrors build_patient_params logic, silent)
params_base = params_default;
params_base.HR_bpm    = clinical.HR_bpm;
params_base.T_cardiac = 60 / clinical.HR_bpm;
params_base.Ts1       = 0.3  * sqrt(params_base.T_cardiac);
params_base.Ts2       = 0.45 * sqrt(params_base.T_cardiac);
params_base.R_systemic = clinical.P_ao_mean_mmHg / clinical.CO_mLs;

P_lv_target    = clinical.P_ao_mean_mmHg * 1.30;
params_base.Emax_lv = P_lv_target / clinical.SV_mL;
params_base.Emin_lv = params_base.Emax_lv * 0.05;
params_base.Emax_rv = params_base.Emax_lv * 0.5;
params_base.Emin_rv = params_base.Emin_lv;

D_pda_m   = clinical.D_shunt_pda_mm * uc.mm_to_m;
A_pda_m2  = pi * (D_pda_m / 2)^2;
Q_pda_est = max(A_pda_m2 * clinical.v_pda_ms * uc.m3s_to_mLs, 0.5);
params_base.R_shunt_pda = max(0.01, min(50, clinical.dP_pda_mmHg / Q_pda_est));
params_base.P_pa_target_mmHg = clinical.P_pa_est_mmHg;

params_base.X0(params_base.idx.P_ao)  = clinical.P_ao_mean_mmHg;
params_base.X0(params_base.idx.P_sys) = clinical.P_ao_mean_mmHg;
params_base.X0(params_base.idx.P_pa)  = clinical.P_pa_est_mmHg;
params_base.X0(params_base.idx.P_pv)  = max(clinical.P_pa_est_mmHg - 5, 3);
params_base.X0(params_base.idx.P_la)  = max(clinical.P_pa_est_mmHg - 7, 3);

fprintf('  R_systemic (baseline): %.4f mmHg·s/mL\n', params_base.R_systemic);
fprintf('  Emax_lv    (baseline): %.4f mmHg/mL\n',   params_base.Emax_lv);
fprintf('  C_sys      (baseline): %.4f mL/mmHg\n\n', params_base.C_sys);

% =========================================================================
%  STEP 3 — BUILD OPTIMIZATION CONFIG STRUCT
% =========================================================================
fprintf('STEP 3: Configuring optimization...\n');

% Extract bounds into separate lb / ub vectors for fmincon
lb = opt_bounds(:, 1);
ub = opt_bounds(:, 2);

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
for k = 1:length(opt_param_names)
    pname = opt_param_names{k};
    switch pname
        case 'R_systemic',    x0(k) = params_base.R_systemic;
        case 'C_sys',         x0(k) = params_base.C_sys;
        case 'Emax_lv',       x0(k) = params_base.Emax_lv;
        case 'stenosis_pct',  x0(k) = default_stenosis_pct;
        case 'coa_length_mm', x0(k) = default_coa_length_mm;
        case 'R_shunt_pda',   x0(k) = params_base.R_shunt_pda;
        case 'C_ao',          x0(k) = params_base.C_ao;
        case 'R_pa',          x0(k) = params_base.R_pa;
        otherwise,            x0(k) = (lb(k) + ub(k)) / 2;
    end
    % Clamp x0 to bounds
    x0(k) = max(lb(k), min(ub(k), x0(k)));
end

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
OPT_HISTORY_X   = x0';         % 1 × D
OPT_HISTORY_J   = J_baseline;  % 1 × 1
OPT_ITER_COUNT  = 0;

obj_func = @(x) objective_lbfgsb_pda_coa(x, opt_config);

t_opt_start = tic;

[x_opt, J_opt, exitflag, output_struct] = fmincon( ...
    obj_func, x0, ...      % objective and initial point
    [], [], [], [], ...    % no linear constraints
    lb, ub, ...            % bounds
    [], ...                % no nonlinear constraints
    fmincon_opts);

t_opt_elapsed = toc(t_opt_start);

fprintf('\n  Optimization finished in %.1f s (%d iterations, %d func evals)\n', ...
    t_opt_elapsed, output_struct.iterations, output_struct.funcCount);
fprintf('  Exit flag: %d  (%s)\n\n', exitflag, exit_flag_message(exitflag));

% Retrieve history recorded by callback
x_history = OPT_HISTORY_X;   % n_iter × D
J_history  = OPT_HISTORY_J;  % n_iter × 1

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
mat_path = fullfile(results_dir, 'optimization_workspace.mat');
save(mat_path, 'x_opt', 'x0', 'J_opt', 'J_baseline', 'J_final', ...
    'x_history', 'J_history', 'opt_config', 'clinical', 'params_opt', ...
    'stenosis_opt', 'coa_length_opt', 'baseline_outputs', 'opt_outputs');
fprintf('  Saved: %s\n\n', mat_path);

% =========================================================================
%  STEP 8 — GENERATE PLOTS
% =========================================================================
fprintf('STEP 8: Generating plots...\n');

plot_optimization_results(baseline_outputs, opt_outputs, clinical, ...
    x_history, J_history, opt_param_names, results_dir);

fprintf('\n');

% =========================================================================
%  STEP 9 — FINAL CLINICAL REPORT
% =========================================================================
fprintf('=========================================================================\n');
fprintf('   OPTIMIZATION RESULTS: CoA CLINICAL OUTPUTS\n');
fprintf('   Patient: %s\n', clinical.patient_id);
fprintf('=========================================================================\n\n');

fprintf('  %-30s  %12s  %12s  %12s\n', 'Target', 'Clinical', 'Pre-Opt', 'Post-Opt');
fprintf('  %s\n', repmat('-', 1, 70));

% Helper for safe display
def_nan = @(s, f) get_val_safe(s, f);

fprintf('  %-30s  %12.1f  %12.1f  %12.1f\n', 'MAP (mmHg)', ...
    clinical.P_ao_mean_mmHg, def_nan(baseline_outputs,'P_ao_mean'), def_nan(opt_outputs,'P_ao_mean'));
fprintf('  %-30s  %12.1f  %12.1f  %12.1f\n', 'SBP (mmHg)', ...
    clinical.P_ao_sys_mmHg, def_nan(baseline_outputs,'P_ao_sys'), def_nan(opt_outputs,'P_ao_sys'));
fprintf('  %-30s  %12.1f  %12.1f  %12.1f\n', 'DBP (mmHg)', ...
    clinical.P_ao_dia_mmHg, def_nan(baseline_outputs,'P_ao_dia'), def_nan(opt_outputs,'P_ao_dia'));
fprintf('  %-30s  %12.2f  %12.2f  %12.2f\n', 'SV (mL)', ...
    clinical.SV_mL, def_nan(baseline_outputs,'SV_lv'), def_nan(opt_outputs,'SV_lv'));
fprintf('\n');
fprintf('  %-30s  %12.1f  %12.1f  %12.1f\n', 'dP_CoA_peak (mmHg)', ...
    clinical.dP_coa_mmHg, def_nan(baseline_outputs,'DeltaP_coa_peak'), def_nan(opt_outputs,'DeltaP_coa_peak'));
fprintf('  %-30s  %12s  %12.1f  %12.1f\n', 'dP_CoA_mean_sys (mmHg)', ...
    '(Doppler)', def_nan(baseline_outputs,'DeltaP_coa_mean_sys'), def_nan(opt_outputs,'DeltaP_coa_mean_sys'));
fprintf('  %-30s  %12s  %12.3f  %12.3f\n', 'Q_CoA / Q_total', ...
    '—', def_nan(baseline_outputs,'Q_coa_fraction'), def_nan(opt_outputs,'Q_coa_fraction'));
fprintf('\n');
fprintf('  %-30s  %12s  %12s  %12s\n', 'Predicted Severity', ...
    '(Echo)', ...
    upper(get_str_safe(baseline_outputs,'predicted_CoA_severity')), ...
    upper(get_str_safe(opt_outputs,'predicted_CoA_severity')));
fprintf('\n');
fprintf('  Objective J:    Baseline = %.6f  |  Final = %.6f\n', J_baseline, J_final);
fprintf('  Improvement:    %.2f%%\n', 100*(J_baseline - J_final)/max(J_baseline,eps));
fprintf('\n');
fprintf('  Optimized CoA geometry:\n');
fprintf('    stenosis_pct  = %.2f%%\n', stenosis_opt);
fprintf('    coa_length_mm = %.2f mm\n\n', coa_length_opt);
fprintf('  Results saved to: %s/\n', results_dir);
fprintf('=========================================================================\n');


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function stop = optimization_output_callback(x, optimValues, state)
% Records x and J at each fmincon iteration into global history arrays.
    global OPT_HISTORY_X OPT_HISTORY_J OPT_ITER_COUNT;
    stop = false;
    if strcmp(state, 'iter')
        OPT_ITER_COUNT  = OPT_ITER_COUNT + 1;
        OPT_HISTORY_J   = [OPT_HISTORY_J;   optimValues.fval];
        OPT_HISTORY_X   = [OPT_HISTORY_X;   x(:)'];
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
