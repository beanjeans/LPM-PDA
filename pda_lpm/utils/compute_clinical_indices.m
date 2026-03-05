function indices = compute_clinical_indices(t_sol, X_sol, params, clinical, scenario_label)
% COMPUTE_CLINICAL_INDICES
% -----------------------------------------------------------------------
% Extracts and computes all derived haemodynamic indices from the ODE
% solution at periodic steady state. Separates model outputs (model.*)
% from clinical measurements (clinical.*) per Guardrail §3.11.
%
% INPUTS:
%   t_sol          - time vector (steady-state reporting cycles)    [s]
%   X_sol          - state matrix (n_points × n_states)
%   params         - parameter struct (includes idx, scenario_coa)
%   clinical       - clinical measurement struct (from load_patient_data)
%   scenario_label - string label for display (e.g., 'PDA only', 'CoA 75%')
%
% OUTPUTS:
%   indices  - struct with fields:
%       .model.*     — model-derived quantities
%       .clinical.*  — reference clinical values (copied, not overwritten)
%       .validation.*— comparison metrics
%
% ASSUMPTIONS:
%   - Last cardiac cycle extracted for all waveform analyses
%   - LV volume reconstructed from P_lv and E_lv: V_lv = P_lv/E_lv + V0_lv
%   - Qp/Qs computed from pulmonary vs systemic flow integrals
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%       Eqs. (21)–(23): MAP, SV, CO.
%   [2] Keshavarz-Motamed et al. (2011). J Biomech 44:2817–2825.
%       CoA gradient metrics.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

idx    = params.idx;
mLs_to_Lmin = 60 / 1000;  % [L/min per mL/s] — named conversion per Guardrail §6.3

%% Extract last complete cardiac cycle for waveform indices
T      = params.T_cardiac;                           % [s]
t_last = t_sol(end) - T;
cycle_mask = t_sol >= t_last;
t_cyc  = t_sol(cycle_mask);
X_cyc  = X_sol(cycle_mask, :);

%% -----------------------------------------------------------------------
%  MODEL PRESSURES
% -----------------------------------------------------------------------
model.P_ra_mean   = mean(X_cyc(:, idx.P_ra));          % [mmHg]
model.P_rv_sys    = max(X_cyc(:, idx.P_rv));            % [mmHg]
model.P_rv_dia    = min(X_cyc(:, idx.P_rv));            % [mmHg]
model.P_pa_sys    = max(X_cyc(:, idx.P_pa));            % [mmHg]
model.P_pa_dia    = min(X_cyc(:, idx.P_pa));            % [mmHg]
model.P_pa_mean   = mean(X_cyc(:, idx.P_pa));           % [mmHg]
model.P_la_mean   = mean(X_cyc(:, idx.P_la));           % [mmHg]
model.P_lv_sys    = max(X_cyc(:, idx.P_lv));            % [mmHg]
model.P_lv_dia    = min(X_cyc(:, idx.P_lv));            % [mmHg]
model.P_ao_sys    = max(X_cyc(:, idx.P_ao));            % [mmHg]
model.P_ao_dia    = min(X_cyc(:, idx.P_ao));            % [mmHg]
model.P_ao_mean   = mean(X_cyc(:, idx.P_ao));           % [mmHg]

%% -----------------------------------------------------------------------
%  CoA-SPECIFIC PRESSURES (only when CoA scenario active)
% -----------------------------------------------------------------------
if params.scenario_coa && isfield(idx, 'P_ao_dist')
    model.P_ao_dist_sys  = max(X_cyc(:, idx.P_ao_dist)); % [mmHg]
    model.P_ao_dist_dia  = min(X_cyc(:, idx.P_ao_dist)); % [mmHg]
    model.P_ao_dist_mean = mean(X_cyc(:, idx.P_ao_dist));% [mmHg]

    DeltaP_coa_trace     = X_cyc(:, idx.P_ao) - X_cyc(:, idx.P_ao_dist);
    model.DeltaP_coa_mean = mean(DeltaP_coa_trace);       % [mmHg]
    model.DeltaP_coa_peak = max(DeltaP_coa_trace);        % [mmHg]

    % CoA severity classification (clinical standard — Ref [2])
    if model.DeltaP_coa_mean < 10
        model.CoA_severity_gradient = 'Mild (<10 mmHg)';
    elseif model.DeltaP_coa_mean < 20
        model.CoA_severity_gradient = 'Moderate (10–20 mmHg)';
    else
        model.CoA_severity_gradient = 'Severe (>20 mmHg)';
    end

    % CoA flow
    model.Q_coa_mean  = mean(X_cyc(:, idx.Q_coa));       % [mL/s]
    model.Q_coa_peak  = max(X_cyc(:, idx.Q_coa));        % [mL/s]
else
    model.DeltaP_coa_mean        = 0;
    model.DeltaP_coa_peak        = 0;
    model.CoA_severity_gradient  = 'N/A (PDA-only scenario)';
    model.Q_coa_mean             = 0;
end

%% -----------------------------------------------------------------------
%  LV VOLUME AND STROKE WORK
% -----------------------------------------------------------------------
% Reconstruct En for last cycle to get E_lv trace
tn_cyc = mod(t_cyc, T);
En_cyc = zeros(size(tn_cyc));
for k = 1:length(tn_cyc)
    if tn_cyc(k) <= params.Ts1
        En_cyc(k) = 1.55 * (tn_cyc(k) / params.Ts1)^2;
    elseif tn_cyc(k) <= params.Ts2
        alpha     = (pi/2) * (tn_cyc(k) - params.Ts1) / (params.Ts2 - params.Ts1);
        En_cyc(k) = 1.55 * cos(alpha)^2;
    end
end
E_lv_cyc  = (params.Emax_lv - params.Emin_lv) * En_cyc + params.Emin_lv; % [mmHg/mL]
V_lv_cyc  = X_cyc(:, idx.P_lv) ./ E_lv_cyc + params.V0_lv;               % [mL]

model.V_lv_ed  = max(V_lv_cyc);     % [mL] — end-diastolic volume
model.V_lv_es  = min(V_lv_cyc);     % [mL] — end-systolic volume
model.SV_lv    = model.V_lv_ed - model.V_lv_es;  % [mL] — stroke volume
model.EF_lv    = model.SV_lv / model.V_lv_ed;    % [fraction]

% LV stroke work = ∫ P dV over one cycle (area of P-V loop) [mmHg·mL]
model.SW_lv_mmHg_mL = abs(trapz(V_lv_cyc, X_cyc(:, idx.P_lv)));  % [mmHg·mL]
uc = unit_conversion();
model.SW_lv_J       = model.SW_lv_mmHg_mL * uc.mmHg_mL_to_J;     % [J]

%% -----------------------------------------------------------------------
%  CARDIAC OUTPUT AND Qp/Qs
% -----------------------------------------------------------------------
model.CO_Lmin = model.SV_lv * params.HR_bpm / 1000;   % [L/min]

% Pulmonary flow: integrate Q_pa_pul over one cycle
Q_pa_cyc      = X_cyc(:, idx.Q_pa_pul);               % [mL/s]
SV_pulmonary  = trapz(t_cyc, max(0, Q_pa_cyc));       % [mL/beat]
CO_pul_Lmin   = SV_pulmonary * params.HR_bpm / 1000;  % [L/min]

model.CO_pul_Lmin = CO_pul_Lmin;                       % [L/min]

% Qp/Qs ratio (should be > 1 in L→R PDA shunt)
if model.CO_Lmin > 0
    model.Qp_Qs = model.CO_pul_Lmin / model.CO_Lmin;  % [dimensionless]
else
    model.Qp_Qs = NaN;
end

% PDA shunt flow mean
model.Q_shunt_pda_mean = mean(X_cyc(:, idx.Q_shunt_pda));  % [mL/s]

%% -----------------------------------------------------------------------
%  COPY CLINICAL REFERENCE VALUES (no modification — Guardrail §3.11)
% -----------------------------------------------------------------------
indices.clinical = clinical;

%% -----------------------------------------------------------------------
%  VALIDATION — MODEL vs. CLINICAL
% -----------------------------------------------------------------------
validation.MAP_error_mmHg = model.P_ao_mean - clinical.P_ao_mean_mmHg;
validation.MAP_error_pct  = 100 * abs(validation.MAP_error_mmHg) / clinical.P_ao_mean_mmHg;

validation.SV_error_mL    = model.SV_lv - clinical.SV_mL;
validation.SV_error_pct   = 100 * abs(validation.SV_error_mL) / clinical.SV_mL;

%% -----------------------------------------------------------------------
%  DISPLAY SUMMARY
% -----------------------------------------------------------------------
fprintf('\n=== HAEMODYNAMIC INDICES: %s ===\n', scenario_label);
fprintf('  LV: Sys/Dia = %.1f / %.1f mmHg  |  EF = %.1f%%\n', ...
    model.P_lv_sys, model.P_lv_dia, model.EF_lv * 100);
fprintf('  Ao: Sys/Dia/Mean = %.1f / %.1f / %.1f mmHg\n', ...
    model.P_ao_sys, model.P_ao_dia, model.P_ao_mean);
fprintf('  PA: Sys/Dia/Mean = %.1f / %.1f / %.1f mmHg\n', ...
    model.P_pa_sys, model.P_pa_dia, model.P_pa_mean);
fprintf('  SV: %.2f mL (clinical: %.2f mL)  |  CO: %.2f L/min\n', ...
    model.SV_lv, clinical.SV_mL, model.CO_Lmin);
fprintf('  Qp/Qs = %.2f  |  Q_PDA_mean = %.2f mL/s\n', ...
    model.Qp_Qs, model.Q_shunt_pda_mean);
fprintf('  LV Stroke Work: %.2f mmHg·mL = %.4f J\n', ...
    model.SW_lv_mmHg_mL, model.SW_lv_J);
if params.scenario_coa
    fprintf('  CoA: Mean dP = %.1f mmHg  |  Peak dP = %.1f mmHg  →  %s\n', ...
        model.DeltaP_coa_mean, model.DeltaP_coa_peak, model.CoA_severity_gradient);
end
fprintf('  Validation: MAP error = %.1f mmHg (%.1f%%)  |  SV error = %.1f mL (%.1f%%)\n\n', ...
    validation.MAP_error_mmHg, validation.MAP_error_pct, ...
    validation.SV_error_mL, validation.SV_error_pct);

%% Pack outputs
indices.model      = model;
indices.validation = validation;
indices.scenario   = scenario_label;

end
