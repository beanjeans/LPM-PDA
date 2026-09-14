function [J, sim_outputs] = objective_lbfgsb_pda_coa(x_opt, opt_config)
% OBJECTIVE_LBFGSB_PDA_COA
% -----------------------------------------------------------------------
% Objective function for bounded quasi-Newton (L-BFGS-B-style) optimization
% of the PDA-CoA Lumped Parameter Model — Patient 1 (PDA-only) configuration.
%
% Patient 1 is a PDA-only patient.  No clinical CoA gradient is measured.
% The objective is restricted to MAP and SV only:
%
%   J = w_MAP * ((sim_MAP - clin_MAP) / max(|clin_MAP|, 1.0))^2
%     + w_SV  * ((sim_SV  - clin_SV)  / max(|clin_SV|,  0.1))^2
%     + plausibility_penalties
%
% Active parameters: Emax_lv, stenosis_pct, R_systemic
% Clinical targets:  MAP, SV
% NOT used in objective: SBP, DBP, PP, dP_PDA, dP_CoA, Q_CoA/Q_total
%
% INPUTS:
%   x_opt       - (D_opt × 1) vector of current parameter values,
%                 in the order defined by opt_config.param_names
%   opt_config  - struct with all fixed optimization settings:
%       .param_names    — cell array of optimized parameter names
%       .fixed_params   — params struct with non-optimized fields fixed
%       .clinical       — clinical measurement struct
%       .stenosis_pct   — CoA stenosis [%] (optimized if in param_names)
%       .coa_length_mm  — CoA length   [mm] (optimized if in param_names)
%       .weights        — struct of objective weights
%       .n_warmup       — ODE warm-up cycles
%       .n_report       — ODE reporting cycles
%       .penalty        — large value returned on solver failure
%
% OUTPUTS:
%   J           - scalar objective value (≥ 0)
%   sim_outputs - struct of simulated quantities (empty on failure);
%                 includes objective_breakdown sub-struct
%
% AUTHOR:   Optimization Extension — Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  2.0  — restricted to MAP+SV objective for PDA-only patient
% -----------------------------------------------------------------------

%% 1. Unpack current parameter vector into the params struct
params = opt_config.fixed_params;
clinical = opt_config.clinical;
w = opt_config.weights;

stenosis_pct  = opt_config.stenosis_pct;
coa_length_mm = opt_config.coa_length_mm;

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
            if ~ismember('Emin_lv', opt_config.param_names)
                params.Emin_lv = val * 0.05;   % Fallback 5% ratio if Emin_lv is not independently optimized
            end
            params.Emax_rv = val * 0.5;
            params.Emin_rv = params.Emin_lv;
        case 'Emin_lv'
            params.Emin_lv = val;              % Direct assignment for LV diastolic elastance
            params.Emin_rv = val;
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
    end
end

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
% Patient 1 is PDA-only.  Only MAP and SV are direct clinical targets.
% Normalization uses max(|clin|, floor) to avoid division by tiny values.
% A term is only added when: weight > 0, clinical value is finite and
% non-zero, and the simulated value is finite.
%
% SBP, DBP, PP, dP_PDA, dP_CoA are NOT included (weights = 0).
% CoA gradient (dP_coa_mmHg = 0 in CSV) is not a measured clinical value
% for this patient and must not be used as an optimization target.
J = 0;
objective_breakdown = struct();

% --- Mean Arterial Pressure (MAP) ---
clin_MAP = clinical.P_ao_mean_mmHg;
sim_MAP  = m.P_ao_mean;
if w.MAP > 0 && isfinite(clin_MAP) && clin_MAP ~= 0 && isfinite(sim_MAP)
    norm_MAP = max(abs(clin_MAP), 1.0);
    err_MAP  = (sim_MAP - clin_MAP) / norm_MAP;
    contrib_MAP = w.MAP * err_MAP^2;
    J = J + contrib_MAP;
    objective_breakdown.MAP.error = sim_MAP - clin_MAP;
    objective_breakdown.MAP.weighted_contribution = contrib_MAP;
end

% --- Stroke Volume (SV) ---
clin_SV = clinical.SV_mL;
sim_SV  = m.SV_lv;
if w.SV > 0 && isfinite(clin_SV) && clin_SV ~= 0 && isfinite(sim_SV)
    norm_SV = max(abs(clin_SV), 0.1);
    err_SV  = (sim_SV - clin_SV) / norm_SV;
    contrib_SV = w.SV * err_SV^2;
    J = J + contrib_SV;
    objective_breakdown.SV.error = sim_SV - clin_SV;
    objective_breakdown.SV.weighted_contribution = contrib_SV;
end

% --- Systolic Blood Pressure (SBP) ---
clin_SBP = clinical.P_ao_sys_mmHg;
sim_SBP  = m.P_ao_sys;
if w.SBP > 0 && isfinite(clin_SBP) && clin_SBP ~= 0 && isfinite(sim_SBP)
    norm_SBP = max(abs(clin_SBP), 1.0);
    err_SBP  = (sim_SBP - clin_SBP) / norm_SBP;
    contrib_SBP = w.SBP * err_SBP^2;
    J = J + contrib_SBP;
    objective_breakdown.SBP.error = sim_SBP - clin_SBP;
    objective_breakdown.SBP.weighted_contribution = contrib_SBP;
end

% --- Diastolic Blood Pressure (DBP) ---
clin_DBP = clinical.P_ao_dia_mmHg;
sim_DBP  = m.P_ao_dia;
if w.DBP > 0 && isfinite(clin_DBP) && clin_DBP ~= 0 && isfinite(sim_DBP)
    norm_DBP = max(abs(clin_DBP), 1.0);
    err_DBP  = (sim_DBP - clin_DBP) / norm_DBP;
    contrib_DBP = w.DBP * err_DBP^2;
    J = J + contrib_DBP;
    objective_breakdown.DBP.error = sim_DBP - clin_DBP;
    objective_breakdown.DBP.weighted_contribution = contrib_DBP;
end

% --- Pulse Pressure (PP = SBP - DBP) ---
clin_PP = clin_SBP - clin_DBP;
sim_PP  = sim_SBP - sim_DBP;
if w.PP > 0 && isfinite(clin_PP) && clin_PP ~= 0 && isfinite(sim_PP)
    norm_PP = max(abs(clin_PP), 1.0);
    err_PP  = (sim_PP - clin_PP) / norm_PP;
    contrib_PP = w.PP * err_PP^2;
    J = J + contrib_PP;
    objective_breakdown.PP.error = sim_PP - clin_PP;
    objective_breakdown.PP.weighted_contribution = contrib_PP;
end

% --- Trans-PDA Pressure Gradient (dP_PDA) ---
clin_dP_PDA = clinical.dP_pda_mmHg;
sim_dP_PDA  = m.P_ao_mean - m.P_pa_mean;
if w.dP_PDA > 0 && isfinite(clin_dP_PDA) && clin_dP_PDA ~= 0 && isfinite(sim_dP_PDA)
    norm_dP_PDA = max(abs(clin_dP_PDA), 0.5);
    err_dP_PDA  = (sim_dP_PDA - clin_dP_PDA) / norm_dP_PDA;
    contrib_dP_PDA = w.dP_PDA * err_dP_PDA^2;
    J = J + contrib_dP_PDA;
    objective_breakdown.dP_PDA.error = sim_dP_PDA - clin_dP_PDA;
    objective_breakdown.dP_PDA.weighted_contribution = contrib_dP_PDA;
end

% --- CoA Peak Gradient (dP_CoA_peak) ---
clin_dP_CoA = clinical.dP_coa_mmHg;
sim_dP_CoA  = m.DeltaP_coa_peak;
if w.dP_CoA_peak > 0 && isfinite(clin_dP_CoA) && clin_dP_CoA > 0 && isfinite(sim_dP_CoA)
    norm_dP_CoA = max(abs(clin_dP_CoA), 1.0);
    err_dP_CoA  = (sim_dP_CoA - clin_dP_CoA) / norm_dP_CoA;
    contrib_dP_CoA = w.dP_CoA_peak * err_dP_CoA^2;
    J = J + contrib_dP_CoA;
    objective_breakdown.dP_CoA_peak.error = sim_dP_CoA - clin_dP_CoA;
    objective_breakdown.dP_CoA_peak.weighted_contribution = contrib_dP_CoA;
end

objective_breakdown.total = J;

%% 7. Pack simulated outputs for caller inspection
% -----------------------------------------------------------------------
% Primary optimization outputs
sim_outputs.P_ao_mean = m.P_ao_mean;
sim_outputs.SV_lv     = m.SV_lv;

% Additional simulated values (not optimization targets for Patient 1)
sim_outputs.P_ao_sys             = m.P_ao_sys;
sim_outputs.P_ao_dia             = m.P_ao_dia;
sim_outputs.PP                   = m.P_ao_sys - m.P_ao_dia;
sim_outputs.P_pa_mean            = m.P_pa_mean;
sim_outputs.CO_Lmin              = m.CO_Lmin;
sim_outputs.EF_lv                = m.EF_lv;
sim_outputs.DeltaP_coa_peak      = m.DeltaP_coa_peak;
sim_outputs.DeltaP_coa_mean_sys  = m.DeltaP_coa_mean_sys;
sim_outputs.Q_coa_fraction       = m.Q_coa_fraction;
sim_outputs.predicted_CoA_severity = m.predicted_CoA_severity;
sim_outputs.pda_modifier_note    = m.pda_modifier_note;

% Objective metadata
sim_outputs.J                    = J;
sim_outputs.objective_breakdown  = objective_breakdown;
sim_outputs.t_sol  = t_sol;
sim_outputs.X_sol  = X_sol;
sim_outputs.params = params_coa;

end
