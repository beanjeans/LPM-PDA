function params = build_patient_params(clinical, params_default)
% BUILD_PATIENT_PARAMS
% -----------------------------------------------------------------------
% Calibrates the LPM parameter struct to a specific PDA patient by
% overriding default values with patient-derived quantities.
%
% Calibration strategy:
%   1. Total systemic resistance R_systemic = MAP / CO
%   2. PDA shunt resistance R_shunt_pda is derived from clinical PDA
%      diameter and Doppler pressure gradient (Bernoulli)
%   3. Ventricular elastance (Emax_lv) scaled to match measured MAP and SV
%   4. Cardiac timing scaled to patient HR
%
% INPUTS:
%   clinical       - clinical struct from load_patient_data.m
%   params_default - default parameter struct from default_parameters.m
%
% OUTPUTS:
%   params  - patient-specific parameter struct ready for integration
%
% ASSUMPTIONS:
%   - Blood is Newtonian and incompressible (μ = 0.004 Pa·s, ρ = 1060 kg/m³)
%   - PDA resistance estimated from Poiseuille + Bernoulli; L_shunt_pda fixed
%   - Compliance values retain default neonate scaling (no patient-specific
%     echo-derived compliance available in the dataset)
%
% SIGN CONVENTIONS:
%   - Q_shunt_pda > 0 : left-to-right (aorta → PA), physiological in PDA
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%   [2] Keshavarz-Motamed et al. (2011). J Biomech 44:2817–2825. (CoA physics)
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

uc     = unit_conversion();
params = params_default;   % Start from validated defaults

%% 1. Cardiac Timing (patient HR)
params.HR_bpm    = clinical.HR_bpm;               % [bpm]
params.T_cardiac = 60 / clinical.HR_bpm;          % [s]
params.Ts1       = 0.3  * sqrt(params.T_cardiac); % [s] — Ref [1] Eq. (7)
params.Ts2       = 0.45 * sqrt(params.T_cardiac); % [s]

%% 2. Systemic Vascular Resistance
% R_systemic = MAP / CO_mLs   [mmHg·s/mL]
% Source: Ohm's law analogy (Ref [1] Eq. 3)
params.R_systemic = clinical.P_ao_mean_mmHg / clinical.CO_mLs;  % [mmHg·s/mL]

%% 3. Ventricular Elastance Calibration
% Target peak LV systolic pressure ≈ 1.30 × MAP (overcomes systemic resistance
% and provides driving pressure for SV). Factor 1.30 from Ref [1] clinical data.
P_lv_sys_target_mmHg = clinical.P_ao_mean_mmHg * 1.30;   % [mmHg]
SV_mL                = clinical.SV_mL;                    % [mL]

% Emax_lv = P_lv_sys_target / SV_mL   [mmHg/mL]
% Derived from E(t) = P(t) / (V(t) - V0): at end-systole, P = Emax*(SV)
params.Emax_lv  = P_lv_sys_target_mmHg / SV_mL;          % [mmHg/mL]
params.Emin_lv  = params.Emax_lv * 0.05;                  % [mmHg/mL] — 5% of Emax, Ref [1]

% RV elastance: ~50% of LV in neonates (elevated PVR at birth) — Ref [1]
params.Emax_rv  = params.Emax_lv * 0.5;                   % [mmHg/mL]
params.Emin_rv  = params.Emin_lv;                          % [mmHg/mL]

%% 4. PDA Shunt Resistance
% Method: use Doppler-derived pressure gradient and estimated flow
% R_shunt_pda ≈ dP_pda / Q_pda_est
%
% Q_pda estimated from Bernoulli continuity:
%   Q_pda ≈ A_pda × v_pda    [m³/s] → [mL/s]
D_pda_m   = clinical.D_shunt_pda_mm * uc.mm_to_m;        % [mm] → [m]
A_pda_m2  = pi * (D_pda_m / 2)^2;                        % [m²]
Q_pda_est_m3s = A_pda_m2 * clinical.v_pda_ms;            % [m³/s]
Q_pda_est_mLs = Q_pda_est_m3s * uc.m3s_to_mLs;          % [mL/s]

% Guard against zero flow estimate
if Q_pda_est_mLs < 0.01
    Q_pda_est_mLs = 0.5;   % [mL/s] — conservative fallback for tiny PDA
    warning('BUILD_PATIENT_PARAMS: Q_pda_est near zero; using fallback 0.5 mL/s');
end

dP_pda_mmHg = clinical.dP_pda_mmHg;                      % [mmHg]

% R_shunt_pda = dP / Q   [mmHg·s/mL]
params.R_shunt_pda = dP_pda_mmHg / Q_pda_est_mLs;        % [mmHg·s/mL]

% Clamp to physically meaningful range [0.01, 50] mmHg·s/mL
params.R_shunt_pda = max(0.01, min(50, params.R_shunt_pda));

%% 5. PA Pressure Target (for elastance of RV)
% Estimated from clinical data: P_pa ≈ P_ao - dP_pda (L→R shunt)
% Used to set appropriate RV output pressure target
params.P_pa_target_mmHg = clinical.P_pa_est_mmHg;        % [mmHg]

%% 6. Update Initial Conditions for This Patient
% Scale aortic pressure to patient MAP
params.X0(params.idx.P_ao)  = clinical.P_ao_mean_mmHg;   % [mmHg]
params.X0(params.idx.P_sys) = clinical.P_ao_mean_mmHg;   % [mmHg]
params.X0(params.idx.P_pa)  = clinical.P_pa_est_mmHg;    % [mmHg]
params.X0(params.idx.P_pv)  = max(clinical.P_pa_est_mmHg - 5, 3); % [mmHg]
params.X0(params.idx.P_la)  = max(clinical.P_pa_est_mmHg - 7, 3); % [mmHg]

%% 7. Display calibrated key parameters
fprintf('--- CALIBRATED PATIENT PARAMETERS ---\n');
fprintf('  HR: %d bpm  |  T_cardiac: %.3f s\n', params.HR_bpm, params.T_cardiac);
fprintf('  R_systemic:   %.3f mmHg·s/mL\n', params.R_systemic);
fprintf('  Emax_lv:      %.3f mmHg/mL\n', params.Emax_lv);
fprintf('  R_shunt_pda:  %.3f mmHg·s/mL  (Q_pda_est: %.2f mL/s)\n', ...
    params.R_shunt_pda, Q_pda_est_mLs);
fprintf('  P_pa_target:  %.1f mmHg\n\n', params.P_pa_target_mmHg);

end
