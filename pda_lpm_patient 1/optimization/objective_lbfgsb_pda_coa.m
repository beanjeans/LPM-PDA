function [J, sim_outputs] = objective_lbfgsb_pda_coa(x_opt, opt_config)
% OBJECTIVE_LBFGSB_PDA_COA
% -----------------------------------------------------------------------
% Objective function for bounded quasi-Newton (L-BFGS-B-style) optimization
% of the PDA-CoA Lumped Parameter Model.
%
% Evaluates a weighted, normalized least-squares error between the
% model-simulated haemodynamic quantities and the patient's clinical
% measurements:
%
%   J = Σ  w_i × ( (sim_i − clin_i) / clin_i )²
%
% The function is designed to be passed to fmincon (MATLAB's built-in
% bounded quasi-Newton solver, which uses an L-BFGS-B-style algorithm
% when called with bounds and the 'quasi-newton' sub-problem solver).
%
% INPUTS:
%   x_opt       - (D_opt × 1) vector of current parameter values,
%                 in the order defined by opt_config.param_names
%   opt_config  - struct with all fixed optimization settings:
%       .param_names    — cell array of optimized parameter names
%       .param_lb       — lower bounds (same order as param_names)
%       .param_ub       — upper bounds
%       .fixed_params   — params struct with non-optimized fields fixed
%       .clinical       — clinical measurement struct
%       .stenosis_pct   — CoA stenosis [%] (optimized if in param_names)
%       .coa_length_mm  — CoA length   [mm] (optimized if in param_names)
%       .weights        — struct of objective weights (w_MAP, w_SBP, etc.)
%       .n_warmup       — ODE warm-up cycles
%       .n_report       — ODE reporting cycles
%       .penalty        — large value returned on solver failure
%
% OUTPUTS:
%   J           - scalar objective value (≥ 0)
%   sim_outputs - struct of simulated quantities (empty on failure)
%
% REFERENCES:
%   [1] Nocedal J & Wright SJ (2006). Numerical Optimization, 2nd ed.
%       Springer. Chapter 7 (L-BFGS).
%   [2] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%
% AUTHOR:   Optimization Extension — Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

%% 1. Unpack current parameter vector into the params struct
% -----------------------------------------------------------------------
params = opt_config.fixed_params;           % Start from fixed baseline
clinical = opt_config.clinical;
w = opt_config.weights;

% Defaults for CoA geometry (may be overridden below)
stenosis_pct  = opt_config.stenosis_pct;
coa_length_mm = opt_config.coa_length_mm;

% Apply each optimized parameter to the struct
for k = 1:length(opt_config.param_names)
    pname = opt_config.param_names{k};
    val   = x_opt(k);

    switch pname
        case 'R_systemic'
            params.R_systemic = val;
        case 'C_sys'
            params.C_sys = val;
        case 'Emax_lv'
            params.Emax_lv = val;
            params.Emin_lv = val * 0.05;   % Maintain 5% diastolic ratio
            params.Emax_rv = val * 0.5;    % RV ≈ 50% LV (neonatal reference)
            params.Emin_rv = params.Emin_lv;
        case 'stenosis_pct'
            stenosis_pct = val;
        case 'coa_length_mm'
            coa_length_mm = val;
        case 'R_shunt_pda'
            params.R_shunt_pda = val;
        case 'C_ao'
            params.C_ao = val;
        case 'R_pa'
            params.R_pa = val;
        otherwise
            % Unknown parameter name — silently ignore (safety)
    end
end

%% 2. Physical validity guard
% -----------------------------------------------------------------------
% Return heavy penalty if any parameter reaches a non-physical state.
% This prevents the solver from exploring negative resistance/compliance.
if params.R_systemic  <= 0 || params.C_sys      <= 0 || ...
   params.Emax_lv     <= 0 || params.C_ao       <= 0 || ...
   params.R_pa        <= 0 || params.R_shunt_pda <= 0 || ...
   stenosis_pct       <= 0 || stenosis_pct       >= 100 || ...
   coa_length_mm      <= 0
    J          = opt_config.penalty;
    sim_outputs = [];
    return;
end

%% 3. Build CoA parameter struct (calls existing, unmodified function)
% -----------------------------------------------------------------------
try
    % Suppress console output from build_coa_params during iteration
    evalc('params_coa = build_coa_params(params, clinical, stenosis_pct, coa_length_mm);');
catch ME
    % build_coa_params validation failed (e.g., stenosis out of range)
    J          = opt_config.penalty;
    sim_outputs = [];
    return;
end

%% 4. Integrate ODE system (calls existing, unmodified solver)
% -----------------------------------------------------------------------
try
    rhs = @(t, X) system_rhs_pda_coa(t, X, params_coa);

    solver_opts = odeset( ...
        'RelTol',  1e-5, ...
        'AbsTol',  1e-7, ...
        'MaxStep', params_coa.T_cardiac / 50);

    T_cyc  = params_coa.T_cardiac;
    tspan  = [0, (opt_config.n_warmup + opt_config.n_report) * T_cyc];

    [t_full, X_full] = ode15s(rhs, tspan, params_coa.X0, solver_opts);

    % Extract steady-state reporting window
    t_start_report = opt_config.n_warmup * T_cyc;
    mask   = t_full >= t_start_report;
    t_sol  = t_full(mask);
    X_sol  = X_full(mask, :);

    % Safety check: need enough time points for meaningful statistics
    if length(t_sol) < 10
        J          = opt_config.penalty;
        sim_outputs = [];
        return;
    end

catch ME
    % ODE solver diverged / failed
    J          = opt_config.penalty;
    sim_outputs = [];
    return;
end

%% 5. Compute clinical indices (calls existing, unmodified function)
% -----------------------------------------------------------------------
try
    evalc('indices = compute_clinical_indices(t_sol, X_sol, params_coa, clinical, ''OPT'');');
catch ME
    J          = opt_config.penalty;
    sim_outputs = [];
    return;
end

m = indices.model;   % Shorthand to simulated model outputs

%% 6. Compute weighted normalized least-squares objective
% -----------------------------------------------------------------------
% J = Σ  w_i × ( (sim_i − clin_i) / clin_i )²
%
% Each term is normalized by the clinical reference value so that all
% outputs are dimensionless and on a comparable scale regardless of units.
% Weights allow prioritizing clinically important targets.

J = 0;

% --- Mean Arterial Pressure ---
% Clinical reference: clinical.P_ao_mean_mmHg
if w.MAP > 0 && clinical.P_ao_mean_mmHg > 0
    err_MAP = (m.P_ao_mean - clinical.P_ao_mean_mmHg) / clinical.P_ao_mean_mmHg;
    J = J + w.MAP * err_MAP^2;
end

% --- Systolic Blood Pressure ---
if w.SBP > 0 && clinical.P_ao_sys_mmHg > 0
    err_SBP = (m.P_ao_sys - clinical.P_ao_sys_mmHg) / clinical.P_ao_sys_mmHg;
    J = J + w.SBP * err_SBP^2;
end

% --- Diastolic Blood Pressure ---
if w.DBP > 0 && clinical.P_ao_dia_mmHg > 0
    err_DBP = (m.P_ao_dia - clinical.P_ao_dia_mmHg) / clinical.P_ao_dia_mmHg;
    J = J + w.DBP * err_DBP^2;
end

% --- Stroke Volume ---
if w.SV > 0 && clinical.SV_mL > 0
    err_SV = (m.SV_lv - clinical.SV_mL) / clinical.SV_mL;
    J = J + w.SV * err_SV^2;
end

% --- Pulse Pressure (SBP − DBP) --- identifiable via C_ao ---
% PP = SBP − DBP ≈ SV / C_ao. This term is active only when C_ao is being
% optimized. If weights.PP = 0 (C_ao fixed), this block is automatically skipped.
% Clinical target: from measured SSAP and SDAP.
if isfield(w, 'PP') && w.PP > 0
    clin_PP = clinical.P_ao_sys_mmHg - clinical.P_ao_dia_mmHg;
    sim_PP  = m.P_ao_sys - m.P_ao_dia;
    if clin_PP > 0
        err_PP = (sim_PP - clin_PP) / clin_PP;
        J = J + w.PP * err_PP^2;
    end
end

% --- PDA pressure gradient (Doppler-derived reference) ---
% clinical.dP_pda_mmHg is the Bernoulli-derived PDA gradient from echo
% Model: pressure difference across the PDA ≈ P_ao_mean − P_pa_mean
if w.dP_PDA > 0 && clinical.dP_pda_mmHg > 0
    sim_dP_pda = m.P_ao_mean - m.P_pa_mean;
    err_dP_PDA = (sim_dP_pda - clinical.dP_pda_mmHg) / clinical.dP_pda_mmHg;
    J = J + w.dP_PDA * err_dP_PDA^2;
end

% --- CoA pressure gradient (Doppler-derived reference, if available) ---
% clinical.dP_coa_mmHg is the measured echo CoA gradient (may be 0 if not measured)
if w.dP_CoA > 0 && isfield(clinical, 'dP_coa_mmHg') && clinical.dP_coa_mmHg > 0
    err_dP_CoA = (m.DeltaP_coa_peak - clinical.dP_coa_mmHg) / clinical.dP_coa_mmHg;
    J = J + w.dP_CoA * err_dP_CoA^2;
end

%% 7. Pack simulated outputs for caller inspection
% -----------------------------------------------------------------------
sim_outputs.P_ao_mean        = m.P_ao_mean;
sim_outputs.P_ao_sys         = m.P_ao_sys;
sim_outputs.P_ao_dia         = m.P_ao_dia;
sim_outputs.PP               = m.P_ao_sys - m.P_ao_dia;
sim_outputs.P_pa_mean        = m.P_pa_mean;
sim_outputs.SV_lv            = m.SV_lv;
sim_outputs.CO_Lmin          = m.CO_Lmin;
sim_outputs.EF_lv            = m.EF_lv;
sim_outputs.DeltaP_coa_peak       = m.DeltaP_coa_peak;
sim_outputs.DeltaP_coa_mean_sys   = m.DeltaP_coa_mean_sys;
sim_outputs.Q_coa_fraction        = m.Q_coa_fraction;
sim_outputs.predicted_CoA_severity= m.predicted_CoA_severity;
sim_outputs.pda_modifier_note     = m.pda_modifier_note;
sim_outputs.J                     = J;
sim_outputs.t_sol  = t_sol;
sim_outputs.X_sol  = X_sol;
sim_outputs.params = params_coa;

end
