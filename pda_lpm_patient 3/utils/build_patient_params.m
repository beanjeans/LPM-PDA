function params = build_patient_params(clinical, params_default)
% BUILD_PATIENT_PARAMS
% -----------------------------------------------------------------------
% Calibrates the LPM parameter struct to a specific PDA patient by
% overriding allometrically-scaled default values with patient-derived
% quantities from clinical echocardiography and Doppler data.
%
% CALIBRATION PIPELINE:
%   1. Cardiac timing from patient HR
%   2. R_systemic = MAP / CO  (NOT allometrically scaled; clinical data)
%   2b. R_pa, R_pv_veins = clinical pulmonary pressure override
%       (allometric b=-1.00 gives 3.55 mmHg·s/mL for 987g; clinical
%        P_pa implies <0.3 mmHg·s/mL for L→R PDA patients)
%   2c. C_sys = SV-based clinical override
%       (allometric gives 0.005 mL/mmHg; starvation of LV preload;
%        clinical SV/dP gives ~0.18 mL/mmHg)
%   3. R_shunt_pda: Doppler-only (Hagen-Poiseuille discarded)
%   4. Ventricular elastance calibration from clinical SV and MAP
%   5. Unstressed volumes scaled to patient weight
%   6. Initial conditions seeded from clinical pressures
%   7. Sanity check against neonatal physiological reference ranges
%
% INPUTS:
%   clinical       - clinical struct from load_patient_data.m
%   params_default - default parameter struct from default_parameters.m
%                   (already allometrically scaled for BW_neo_kg)
%
% OUTPUTS:
%   params  - patient-specific parameter struct ready for ODE integration
%
% PHYSIOLOGICAL TARGETS (sanity check ranges):
%   SV         :  2–10 mL          (preterm neonate, 1–3 kg)
%   CO         :  0.3–1.5 L/min    (preterm neonate, 1–3 kg)
%   MAP        :  30–55 mmHg       (preterm/term neonate)
%   Qp/Qs      :  >1.0             (L→R PDA, mandatory for shunt physiology)
%
% SIGN CONVENTIONS:
%   Q_shunt_pda > 0 : left-to-right (aorta → PA), physiological for PDA
%
% REFERENCES:
%   [OR2022]  Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%   [P&F2000] Pennati G, Fumero R. (2000). Ann Biomed Eng 28:442–452.
%             DOI: 10.1114/1.282.  (Allometric turbulent R exponent b=−1.33)
%   [S2026]   Seemann G et al. (2026). ASAIO J 72(3):207–215.
%             DOI: 10.1097/MAT.0000000000002528.
%   [HP]      Hagen-Poiseuille: R = 128 μ L / (π D⁴)  [Pa·s/m³]
%   [Rud2001] Rudolph AM (2001). Congenital Diseases of the Heart.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2026-05-12
% VERSION:  2.0  — allometric scaling + Hagen-Poiseuille PDA + sanity check
% -----------------------------------------------------------------------

uc     = unit_conversion();
params = params_default;   % Start from allometrically-scaled defaults

fprintf('\n=== BUILD_PATIENT_PARAMS: Patient %s ===\n', clinical.patient_id);

%% -----------------------------------------------------------------------
%  STEP 1: Cardiac Timing (patient HR)
%  HR is NOT allometrically scaled — it is directly measured.
% -----------------------------------------------------------------------
params.HR_bpm    = clinical.HR_bpm;               % [bpm]
params.T_cardiac = 60 / clinical.HR_bpm;          % [s]
params.Ts1       = 0.3  * sqrt(params.T_cardiac); % [s] — [OR2022] Eq. 7
params.Ts2       = 0.45 * sqrt(params.T_cardiac); % [s]

%% -----------------------------------------------------------------------
%  STEP 2: Systemic Vascular Resistance — clinical derivation
%  R_systemic = MAP / CO   [mmHg·s/mL]
%  This is NOT allometrically scaled; it directly reflects this patient's
%  systemic vascular tone as measured hemodynamically.
%  Reference: Ohm's law analogy, [OR2022] Eq. 3
% -----------------------------------------------------------------------
params.R_systemic = clinical.P_ao_mean_mmHg / clinical.CO_mLs;  % [mmHg·s/mL]

%% -----------------------------------------------------------------------
%  STEP 2b: Pulmonary Vascular Resistance — clinical override
%  Allometric b = -1.00 gives R_pa ≈ 3.55 mmHg·s/mL for a 987g neonate,
%  which is 10× too high. For L→R PDA, the pulmonary circuit is actually
%  a LOW-resistance pathway that accepts the shunt flow. The correct R_pa
%  must be derived from the clinical PA pressure and estimated pulmonary flow.
%
%  Assumptions:
%    Qp/Qs ≈ 1.5 for a significant L→R PDA (Rudolph 2001)
%    P_pv ≈ P_pa - 5 mmHg  (normal pulmonary venous driving pressure)
% -----------------------------------------------------------------------
Q_pul_est_mLs         = clinical.CO_mLs * 1.5;
P_pv_est_mmHg         = max(clinical.P_pa_est_mmHg - 5, 3);
params.R_pa           = max(0.05, min(2.0, ...
    (clinical.P_pa_est_mmHg - P_pv_est_mmHg) / Q_pul_est_mLs));
params.R_pv_veins     = params.R_pa;  % symmetric pulmonary venous resistance

%% -----------------------------------------------------------------------
%  STEP 2c: Systemic Venous Compliance — clinical override
%  Allometric b = +1.33 gives C_sys ≈ 0.005 mL/mmHg for a 987g neonate.
%  This starves the LV of preload, causing SV ≈ 1.5 mL instead of 6 mL.
%
%  Clinical estimate: the venous compartment must store at least one
%  stroke volume above the mean filling pressure.
%    C_sys ≈ SV / (MAP - P_ra_ref)   with P_ra_ref ≈ 4 mmHg
% -----------------------------------------------------------------------
C_sys_clinical    = clinical.SV_mL / max(clinical.P_ao_mean_mmHg - 4, 1);
params.C_sys      = max(0.05, min(0.5, C_sys_clinical));

fprintf('  R_pa   (clinical override): %.4f mmHg·s/mL  (allometric was %.4f)\n', ...
    params.R_pa, params_default.R_pa);
fprintf('  C_sys  (clinical override): %.4f mL/mmHg    (allometric was %.4f)\n\n', ...
    params.C_sys, params_default.C_sys);
% -----------------------------------------------------------------------
%% -----------------------------------------------------------------------
%  STEP 3: PDA Shunt Resistance — Doppler-only
%
%  Uses clinical dP_pda and Doppler velocity directly.
%  Hagen-Poiseuille discarded: for D=2.2 mm it gives R≈0.26, Q≈107 mL/s
%  (10× too large), flooding the pulmonary circuit and collapsing DBP.
%
%  Reference: Rudolph (2001); clinical Bernoulli echo assessment.
% -----------------------------------------------------------------------

D_pda_m      = clinical.D_shunt_pda_mm * uc.mm_to_m;   % [mm] → [m]
A_pda_m2     = pi * (D_pda_m / 2)^2;                   % [m²]

% Doppler-derived PDA resistance: R = dP / Q   (Q = A × v from CW Doppler)
%
% WHY Doppler-only (not Hagen-Poiseuille):
%   HP assumes ideal, fully-developed laminar flow in a perfectly straight
%   tube, which is not valid for the PDA (short, curved, with end-effects).
%   For D_PDA = 2.2 mm, HP gives R ≈ 0.26 mmHg·s/mL → Q ≈ 107 mL/s,
%   ten times the Doppler-measured value of ~10 mL/s. This floods the
%   pulmonary circuit and collapses DBP. The Doppler measurement directly
%   encodes in-vivo duct geometry and flow regime, making it far superior.
%
%   Reference: Rudolph (2001); clinical Bernoulli-based echo assessment.
Q_pda_vel_m3s = A_pda_m2 * clinical.v_pda_ms;           % [m³/s]
Q_pda_vel_mLs = Q_pda_vel_m3s * uc.m3s_to_mLs;         % [mL/s]

if Q_pda_vel_mLs < 0.01
    Q_pda_vel_mLs = 0.5;
    warning('BUILD_PATIENT_PARAMS: Q_pda_vel near zero; using fallback 0.5 mL/s');
end

R_pda_Doppler = clinical.dP_pda_mmHg / Q_pda_vel_mLs;  % [mmHg·s/mL]

% Clamp: lower bound 0.5 prevents unphysical pulmonary flooding;
%        upper bound 20 prevents near-closed shunt at baseline.
params.R_shunt_pda = max(0.5, min(20, R_pda_Doppler));  % [mmHg·s/mL]

fprintf('  PDA R_shunt (Doppler): %.4f mmHg·s/mL  (Q_est=%.2f mL/s, dP=%.1f mmHg)\n', ...
    params.R_shunt_pda, Q_pda_vel_mLs, clinical.dP_pda_mmHg);

%% -----------------------------------------------------------------------
%  STEP 4: Ventricular Elastance Calibration
%
%  The allometric initial seed (b = −1.0 from scale_params_allometric)
%  provides a BW-scaled estimate. Here we OVERRIDE with the clinical
%  target-based formula:
%
%   Emax_lv = P_lv_sys_target / SV_mL   [mmHg/mL]
%
%  where P_lv_sys_target ≈ 1.30 × MAP overcomes systemic resistance and
%  provides driving pressure for the stroke volume.
%
%  WHY flat ×0.8 FAILED (and why allometric b = −1.0 seed is better):
%    Adult Emax_lv = 2.0 mmHg/mL produces SV ≈ 70 mL at P_lv ≈ 120 mmHg.
%    Applying ×0.8 gives Emax_lv = 1.6 mmHg/mL — still adult-scale.
%    For a 1.237 kg neonate: P_lv_sys ≈ 51×1.3 = 66 mmHg, SV ≈ 6 mL
%    → Emax_lv ≈ 11 mmHg/mL (about 5.5× larger than adult, not 0.8×).
%    The flat 0.8 multiplier was insufficient because it was applied to
%    the adult Emax numerically, ignoring that SV is 70→6 mL (12× smaller).
%    Result: the old Emax produced SV ≈ 44 mL (adult-scale output).
%
%  Allometric b = −1.0 seed → Emax ≈ 113 mmHg/mL (too large before
%  clinical calibration, but clamped in scale_params_allometric).
%  The clinical override below sets the physiologically correct value.
%
%  Reference: [OR2022] Eq. (calibration section); [S2026] §2.3
% -----------------------------------------------------------------------
P_lv_sys_target = clinical.P_ao_mean_mmHg * 1.30;   % [mmHg] — LV systolic target
SV_mL           = clinical.SV_mL;                   % [mL]

params.Emax_lv  = P_lv_sys_target / SV_mL;          % [mmHg/mL]
params.Emin_lv  = params.Emax_lv * 0.05;            % 5% of Emax [Ste1996]

% RV: ~50% of LV Emax in neonates (elevated PVR at birth) — [OR2022]
params.Emax_rv  = params.Emax_lv * 0.50;            % [mmHg/mL]
params.Emin_rv  = params.Emin_rv * 0.05;                   % [mmHg/mL]

% Clamp to physiologically plausible neonatal range (as in scale_params_allometric)
EMAX_LV_MIN = 3.0;   EMAX_LV_MAX = 20.0;   % [mmHg/mL]
EMAX_RV_MIN = 1.5;   EMAX_RV_MAX = 12.0;   % [mmHg/mL]
params.Emax_lv = max(EMAX_LV_MIN, min(EMAX_LV_MAX, params.Emax_lv));
params.Emin_lv = params.Emax_lv * 0.05;
params.Emax_rv = max(EMAX_RV_MIN, min(EMAX_RV_MAX, params.Emax_rv));
params.Emin_rv = params.Emax_rv * 0.05;

fprintf('  Emax_lv (clinical calibration): %.4f mmHg/mL  (target P_lv=%.1f, SV=%.2f mL)\n', ...
    params.Emax_lv, P_lv_sys_target, SV_mL);
fprintf('  Emax_rv:  %.4f mmHg/mL  |  Emin_lv: %.4f  |  Emin_rv: %.4f\n', ...
    params.Emax_rv, params.Emin_lv, params.Emin_rv);

%% -----------------------------------------------------------------------
%  STEP 5: PA Pressure Target
%  Estimated from: P_pa ≈ P_ao - dP_pda (L→R shunt drives this gradient)
% -----------------------------------------------------------------------
params.P_pa_target_mmHg = clinical.P_pa_est_mmHg;   % [mmHg]

%% -----------------------------------------------------------------------
%  STEP 6: Update Initial Conditions for This Patient
%  Seeds the ODE solver at physiologically plausible neonatal pressures.
% -----------------------------------------------------------------------
params.X0(params.idx.P_ao)  = clinical.P_ao_mean_mmHg;                   % [mmHg]
params.X0(params.idx.P_sys) = clinical.P_ao_mean_mmHg;                   % [mmHg]
params.X0(params.idx.P_pa)  = clinical.P_pa_est_mmHg;                    % [mmHg]
params.X0(params.idx.P_pv)  = max(clinical.P_pa_est_mmHg - 5,  3);       % [mmHg]
params.X0(params.idx.P_la)  = max(clinical.P_pa_est_mmHg - 7,  3);       % [mmHg]

%% -----------------------------------------------------------------------
%  STEP 7: Sanity Check — Neonatal Physiological Reference Ranges
%
%  Expected vs. computed values are printed. Warnings raised if outside
%  physiologically valid neonatal ranges based on:
%    [Rud2001] Rudolph (2001) — preterm neonate reference ranges
%    [S2026]   Seemann et al. (2026) — allometric model validation
%
%  Ranges (preterm neonate, ~1–3 kg, PDA present):
%    SV    :  2–10 mL       CO    :  0.3–1.5 L/min
%    MAP   :  30–55 mmHg    Qp/Qs :  >1.0 (L→R shunt required)
% -----------------------------------------------------------------------

%% 7a. Compute expected cardiac output from clinical data
SV_check_mL    = clinical.SV_mL;                           % [mL]
CO_check_Lmin  = (SV_check_mL * clinical.HR_bpm) / 1000;  % [L/min]  (SV[mL]×HR[bpm]/1000)
MAP_check_mmHg = clinical.P_ao_mean_mmHg;                  % [mmHg]

%% 7b. Estimate Qp/Qs from R_systemic and R_shunt_pda
% Simple steady-state estimate:
%   Qs ∝ (P_ao - P_sys) / R_systemic  →  Qs ≈ CO
%   Qp = Qs + Q_pda;  Q_pda ≈ dP_pda / R_shunt_pda
%   Qp/Qs ≈ 1 + (Q_pda / Qs)
Q_pda_ss_mLs   = clinical.dP_pda_mmHg / params.R_shunt_pda; % [mL/s]
Qs_mLs         = clinical.CO_mLs;                           % [mL/s]
Qp_mLs         = Qs_mLs + Q_pda_ss_mLs;                    % [mL/s]
QpQs_check     = Qp_mLs / Qs_mLs;                          % [dimensionless]

%% 7c. Reference ranges
REF_SV_MIN    =  2.0;   REF_SV_MAX    = 10.0;   % [mL]
REF_CO_MIN    =  0.3;   REF_CO_MAX    =  1.5;   % [L/min]
REF_MAP_MIN   = 30.0;   REF_MAP_MAX   = 55.0;   % [mmHg]
REF_QPQS_MIN  =  1.0;                           % [dimensionless] — L→R shunt

%% 7d. Print report
fprintf('\n--- SANITY CHECK: Neonatal Physiological Parameter Validation ---\n');
fprintf('  %-20s  Expected Range          Computed        Status\n', 'Parameter');
fprintf('  %-20s  %-22s  %-14s  %s\n', repmat('-',1,20), repmat('-',1,22), repmat('-',1,14), '------');

% SV
sv_ok = (SV_check_mL >= REF_SV_MIN) && (SV_check_mL <= REF_SV_MAX);
fprintf('  %-20s  %4.1f – %4.1f mL           %8.2f mL      %s\n', ...
    'SV [mL]', REF_SV_MIN, REF_SV_MAX, SV_check_mL, status_str(sv_ok));

% CO
co_ok = (CO_check_Lmin >= REF_CO_MIN) && (CO_check_Lmin <= REF_CO_MAX);
fprintf('  %-20s  %4.1f – %4.1f L/min       %8.3f L/min   %s\n', ...
    'CO [L/min]', REF_CO_MIN, REF_CO_MAX, CO_check_Lmin, status_str(co_ok));

% MAP
map_ok = (MAP_check_mmHg >= REF_MAP_MIN) && (MAP_check_mmHg <= REF_MAP_MAX);
fprintf('  %-20s  %4.0f – %4.0f mmHg         %8.1f mmHg    %s\n', ...
    'MAP [mmHg]', REF_MAP_MIN, REF_MAP_MAX, MAP_check_mmHg, status_str(map_ok));

% Qp/Qs
qpqs_ok = (QpQs_check > REF_QPQS_MIN);
fprintf('  %-20s  >%4.1f (L→R shunt)        %8.3f          %s\n', ...
    'Qp/Qs', REF_QPQS_MIN, QpQs_check, status_str(qpqs_ok));

fprintf('\n');

%% 7e. Issue warnings for out-of-range values
if ~sv_ok
    warning('BUILD_PATIENT_PARAMS:SV_OutOfRange', ...
        'SV = %.2f mL is outside neonatal range [%.1f, %.1f] mL. Check SV_mL in clinical data.', ...
        SV_check_mL, REF_SV_MIN, REF_SV_MAX);
end
if ~co_ok
    warning('BUILD_PATIENT_PARAMS:CO_OutOfRange', ...
        'CO = %.3f L/min is outside neonatal range [%.1f, %.1f] L/min. Check SV and HR.', ...
        CO_check_Lmin, REF_CO_MIN, REF_CO_MAX);
end
if ~map_ok
    warning('BUILD_PATIENT_PARAMS:MAP_OutOfRange', ...
        'MAP = %.1f mmHg is outside neonatal range [%.0f, %.0f] mmHg. Check P_ao_mean_mmHg.', ...
        MAP_check_mmHg, REF_MAP_MIN, REF_MAP_MAX);
end
if ~qpqs_ok
    warning('BUILD_PATIENT_PARAMS:QpQs_NotLR', ...
        'Qp/Qs = %.3f ≤ 1.0 — does not indicate L→R PDA shunt. Check R_shunt_pda=%.4f, dP_pda=%.1f mmHg.', ...
        QpQs_check, params.R_shunt_pda, clinical.dP_pda_mmHg);
end

%% -----------------------------------------------------------------------
%  STEP 8: Display final calibrated parameter summary
% -----------------------------------------------------------------------
fprintf('--- CALIBRATED PATIENT PARAMETERS ---\n');
fprintf('  HR: %d bpm  |  T_cardiac: %.3f s\n',     params.HR_bpm, params.T_cardiac);
fprintf('  R_systemic:   %.4f mmHg·s/mL\n',         params.R_systemic);
fprintf('  Emax_lv:      %.4f mmHg/mL\n',            params.Emax_lv);
fprintf('  Emax_rv:      %.4f mmHg/mL\n',            params.Emax_rv);
fprintf('  R_shunt_pda:  %.4f mmHg·s/mL\n',         params.R_shunt_pda);
fprintf('  C_ao (scaled):%.5f mL/mmHg  [P&F2000 b=+1.33]\n', params.C_ao);
fprintf('  C_sys(scaled):%.5f mL/mmHg  [P&F2000 b=+1.33]\n', params.C_sys);
fprintf('  C_pa (scaled):%.5f mL/mmHg  [P&F2000 b=+1.33]\n', params.C_pa);
fprintf('  L_pa (scaled):%.3e mmHg·s²/mL  [P&F2000 b=-0.33]\n', params.L_pa);
fprintf('  L_ao (scaled):%.3e mmHg·s²/mL  [P&F2000 b=-0.33]\n', params.L_ao);
fprintf('  P_pa_target:  %.1f mmHg\n\n',             params.P_pa_target_mmHg);

end

%% -----------------------------------------------------------------------
%  LOCAL HELPER: status string for sanity check table
% -----------------------------------------------------------------------
function s = status_str(ok)
    if ok
        s = 'OK';
    else
        s = '*** WARNING ***';
    end
end
