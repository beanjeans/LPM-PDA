function indices = compute_clinical_indices(t_sol, X_sol, params, clinical, scenario_label)
% COMPUTE_CLINICAL_INDICES
% -----------------------------------------------------------------------
% Extracts and computes all derived haemodynamic indices from the ODE
% solution at periodic steady state. Separates model outputs (model.*)
% from clinical measurements (clinical.*) per Guardrail §3.11.
%
% For CoA scenarios the PRIMARY clinical output is:
%   predicted_CoA_severity  — classified from the simulated ΔP_CoA
%                             (NOT from stenosis_pct directly).
% Thresholds are read from params.severity_thresholds (editable in
% build_coa_params) following ESC guideline consensus:
%     Mild:       mean-systolic ΔP  <  20 mmHg
%     Moderate:   mean-systolic ΔP  20–40 mmHg
%     Severe:     mean-systolic ΔP  >  40 mmHg
%
% PDA parameters are retained and reported as haemodynamic MODIFIERS:
% a large PDA can raise P_ao_dist, reducing observed ΔP_CoA and masking
% true anatomical severity.
%
% INPUTS:
%   t_sol          - time vector (steady-state reporting cycles)    [s]
%   X_sol          - state matrix (n_points × n_states)
%   params         - parameter struct (includes idx, scenario_coa,
%                    severity_thresholds, coa_length_mm, etc.)
%   clinical       - clinical measurement struct (from load_patient_data)
%   scenario_label - string label for display (e.g., 'PDA only', 'CoA 75%')
%
% OUTPUTS:
%   indices  - struct with fields:
%       .model.*     — model-derived quantities (see list below)
%       .clinical.*  — reference clinical values (copied, not overwritten)
%       .validation.*— comparison metrics
%
%   Key model.* fields for CoA scenarios:
%       .stenosis_pct            — anatomical input                [%]
%       .coa_length_mm           — CoA segment length used        [mm]
%       .coa_length_category     — 'discrete/short' or 'long-segment'
%       .DeltaP_coa_peak         — peak pressure gradient         [mmHg]
%       .DeltaP_coa_mean_sys     — mean systolic gradient         [mmHg]
%       .DeltaP_coa_mean         — full-cycle mean gradient       [mmHg]
%       .Q_coa_fraction          — Q_coa_mean / Q_total_mean     [0–1]
%       .P_ao_proximal_mean      — mean proximal aortic pressure  [mmHg]
%       .P_ao_distal_mean        — mean distal aortic pressure    [mmHg]
%       .predicted_CoA_severity  — 'mild'/'moderate'/'severe'
%       .pda_modifier_note       — text note on PDA masking effect
%
% ASSUMPTIONS:
%   - Last cardiac cycle extracted for all waveform analyses
%   - "Systolic" period = time points where P_lv > P_ao (approx. ejection)
%   - LV volume reconstructed from P_lv and E_lv: V_lv = P_lv/E_lv + V0_lv
%   - Qp/Qs computed from pulmonary vs systemic flow integrals
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%       Eqs. (21)–(23): MAP, SV, CO.
%   [2] Keshavarz-Motamed et al. (2011). J Biomech 44:2817–2825.
%       CoA gradient metrics.
%   [3] Baumgartner et al. (2010). Eur Heart J 31(19):2369–2417.
%       ESC guidelines: significant gradient ≥ 20 mmHg; severe > 40 mmHg.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  2.0  — pressure-gradient-based severity classification
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
%  CoA-SPECIFIC OUTPUTS (only when CoA scenario active)
%
%  PRIMARY OUTPUT: predicted_CoA_severity  ← based on ΔP_CoA gradient
%  The anatomical stenosis_pct drives geometry, but severity is REPORTED
%  via the simulated pressure gradient to align with clinical practice.
% -----------------------------------------------------------------------
if params.scenario_coa && isfield(idx, 'P_ao_dist')

    % --- Proximal / distal aortic pressures ----------------------------
    model.P_ao_proximal_sys  = max(X_cyc(:, idx.P_ao));           % [mmHg]
    model.P_ao_proximal_dia  = min(X_cyc(:, idx.P_ao));           % [mmHg]
    model.P_ao_proximal_mean = mean(X_cyc(:, idx.P_ao));          % [mmHg]

    model.P_ao_distal_sys    = max(X_cyc(:, idx.P_ao_dist));      % [mmHg]
    model.P_ao_distal_dia    = min(X_cyc(:, idx.P_ao_dist));      % [mmHg]
    model.P_ao_distal_mean   = mean(X_cyc(:, idx.P_ao_dist));     % [mmHg]

    % Aliases kept for backward compatibility with plot_results
    model.P_ao_dist_sys  = model.P_ao_distal_sys;
    model.P_ao_dist_dia  = model.P_ao_distal_dia;
    model.P_ao_dist_mean = model.P_ao_distal_mean;

    % --- CoA pressure gradient trace ------------------------------------
    DeltaP_coa_trace = X_cyc(:, idx.P_ao) - X_cyc(:, idx.P_ao_dist);  % [mmHg]

    model.DeltaP_coa_peak  = max(DeltaP_coa_trace);    % peak gradient [mmHg]
    model.DeltaP_coa_mean  = mean(DeltaP_coa_trace);   % full-cycle mean [mmHg]

    % --- Mean SYSTOLIC gradient
    %   Approximate systole as frames where P_lv > P_ao (ejection phase)
    sys_mask = X_cyc(:, idx.P_lv) > X_cyc(:, idx.P_ao);
    if any(sys_mask)
        model.DeltaP_coa_mean_sys = mean(DeltaP_coa_trace(sys_mask));  % [mmHg]
    else
        % Fallback: peak-gradient phase (upper quartile of gradient trace)
        thr_q3 = quantile(DeltaP_coa_trace, 0.75);
        model.DeltaP_coa_mean_sys = mean(DeltaP_coa_trace(DeltaP_coa_trace >= thr_q3));
    end

    % --- CoA flow -------------------------------------------------------
    model.Q_coa_mean  = mean(X_cyc(:, idx.Q_coa));     % [mL/s]
    model.Q_coa_peak  = max(X_cyc(:, idx.Q_coa));      % [mL/s]

    % --- CoA / total flow fraction  (Q_coa / Q_total) ------------------
    % Use forward (positive) mean flows to avoid blowup when flows transiently
    % reverse: raw mean can be near-zero or negative, making the ratio explode.
    Q_coa_fwd   = mean(max(0, X_cyc(:, idx.Q_coa)));    % [mL/s]
    Q_ao_fwd    = mean(max(0, X_cyc(:, idx.Q_ao_sys))); % [mL/s]
    Q_total_fwd = Q_coa_fwd + Q_ao_fwd;                 % [mL/s]
    if Q_total_fwd > 1e-6
        model.Q_coa_fraction = Q_coa_fwd / Q_total_fwd;  % guaranteed [0–1]
    else
        model.Q_coa_fraction = NaN;
    end

    % --- Anatomical metadata (pass-through from params) -----------------
    model.stenosis_pct         = params.stenosis_pct;          % [%]
    model.coa_length_mm        = params.coa_length_mm;         % [mm]
    model.coa_length_category  = params.coa_length_category;   % string

    % --- PRIMARY CLINICAL OUTPUT: Severity Classification ---------------
    %   Classification is based on mean SYSTOLIC ΔP_CoA per Ref [3].
    %   Edit params.severity_thresholds in build_coa_params to change.
    if isfield(params, 'severity_thresholds')
        mild_thr = params.severity_thresholds.mild_upper_mmHg;
        mod_thr  = params.severity_thresholds.moderate_upper_mmHg;
    else
        mild_thr = 20;   % mmHg — fallback defaults (Ref [3])
        mod_thr  = 40;   % mmHg
    end

    dp_classify = model.DeltaP_coa_mean_sys;   % gradient used for classification

    if dp_classify < mild_thr
        model.predicted_CoA_severity = 'mild';
    elseif dp_classify < mod_thr
        model.predicted_CoA_severity = 'moderate';
    else
        model.predicted_CoA_severity = 'severe';
    end

    % --- PDA modifier note -----------------------------------------------
    %   PDA size modifies the observed ΔP_CoA. A large PDA raises P_ao_dist
    %   via retrograde collateral flow, reducing the apparent gradient and
    %   potentially causing under-classification of CoA severity.
    if isfield(params, 'R_shunt_pda') && ~isinf(params.R_shunt_pda)
        pda_note = sprintf('PDA present (R_pda=%.3f mmHg·s/mL) — may reduce observed dP_CoA; severity may be underestimated', ...
                           params.R_shunt_pda);
    else
        pda_note = 'PDA absent — ΔP_CoA reflects isolated CoA gradient';
    end
    model.pda_modifier_note = pda_note;

else
    % PDA-only mode — set CoA fields to NaN to prevent misinterpretation as real clinical values.
    % Zero (0) must NOT be used here: a zero gradient is a valid clinical measurement,
    % whereas these patients simply have no CoA module active.
    model.DeltaP_coa_mean        = NaN;
    model.DeltaP_coa_mean_sys    = NaN;
    model.DeltaP_coa_peak        = NaN;
    model.Q_coa_mean             = NaN;
    model.Q_coa_peak             = NaN;
    model.Q_coa_fraction         = NaN;
    model.P_ao_proximal_mean     = model.P_ao_mean;
    model.P_ao_distal_mean       = NaN;
    model.stenosis_pct           = NaN;
    model.coa_length_mm          = NaN;
    model.coa_length_category    = 'N/A (PDA-only)';
    model.predicted_CoA_severity = 'Not applicable';
    model.pda_modifier_note      = 'No CoA simulated — PDA-only mode';
    % Backward-compatible aliases
    model.P_ao_dist_mean         = NaN;
    model.P_ao_dist_sys          = NaN;
    model.P_ao_dist_dia          = NaN;
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
fprintf('  Validation: MAP error = %.1f mmHg (%.1f%%)  |  SV error = %.1f mL (%.1f%%)\n', ...
    validation.MAP_error_mmHg, validation.MAP_error_pct, ...
    validation.SV_error_mL, validation.SV_error_pct);

if params.scenario_coa
    fprintf('\n  --- CoA Clinical Output ---\n');
    fprintf('  Stenosis:          %.0f%%\n',     model.stenosis_pct);
    fprintf('  CoA length:        %.1f mm  [%s]\n', model.coa_length_mm, model.coa_length_category);
    fprintf('  P_ao_proximal:     %.1f mmHg (mean)\n', model.P_ao_proximal_mean);
    fprintf('  P_ao_distal:       %.1f mmHg (mean)\n', model.P_ao_distal_mean);
    fprintf('  ΔP_CoA peak:       %.1f mmHg\n',  model.DeltaP_coa_peak);
    fprintf('  ΔP_CoA mean-sys:   %.1f mmHg\n',  model.DeltaP_coa_mean_sys);
    fprintf('  Q_coa fraction:    %.2f  (Q_coa/Q_total)\n', model.Q_coa_fraction);
    fprintf('  → PREDICTED SEVERITY: %s\n',      upper(model.predicted_CoA_severity));
    fprintf('  PDA modifier: %s\n\n',             model.pda_modifier_note);
end

%% Pack outputs
indices.model      = model;
indices.validation = validation;
indices.scenario   = scenario_label;

end
