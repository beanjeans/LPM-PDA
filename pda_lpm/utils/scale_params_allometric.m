function params_neo = scale_params_allometric(params_adult, BW_neo_kg, BW_adult_kg)
% SCALE_PARAMS_ALLOMETRIC
% -----------------------------------------------------------------------
% Applies body-weight allometric power-law scaling to convert adult-
% reference cardiovascular LPM parameters to neonatal/preterm values.
%
% Scaling law (general form):
%   param_neo = param_adult × (BW_neo / BW_adult) ^ b
%
% Exponents per parameter type are from:
%   [P&F2000] Pennati & Fumero (2000). Ann Biomed Eng 28:442–452.
%             DOI: 10.1114/1.282.
%             Table 1 / Eqs. 4–7 — allometric exponents validated against
%             fetal Doppler data across gestational ages.
%   [S2026]   Seemann et al. (2026). ASAIO J 72(3):207–215.
%             DOI: 10.1097/MAT.0000000000002528.
%             Validates zero-D LPM allometric scaling from birth to adult;
%             Z-score 1.16 for infant simulations.
%
% EXPONENT REFERENCE TABLE:
% ┌──────────────────────────┬──────────┬──────────────────────────────┐
% │ Parameter type           │ Exponent │ Reference                    │
% ├──────────────────────────┼──────────┼──────────────────────────────┤
% │ Viscous resistance (R)   │  −1.00   │ Pennati & Fumero (2000) §3.1 │
% │ Turbulent resistance (R) │  −1.33   │ Pennati & Fumero (2000) §3.2 │
% │ Compliance (C)           │  +1.33   │ Pennati & Fumero (2000) §3.3 │
% │ Inertance (L)            │  −0.33   │ Pennati & Fumero (2000) §3.4 │
% │ Cardiac elastance (E)    │  −1.00   │ Seemann et al. (2026) §2.3   │
% └──────────────────────────┴──────────┴──────────────────────────────┘
%
% PARAMETERS NOT SCALED ALLOMETRICALLY (documented here for traceability):
%   R_systemic  — derived from clinical MAP/CO; allometric scaling would
%                 override patient-specific hemodynamics.  Set in
%                 build_patient_params.m from clinical data.
%   HR          — patient-measured; cardiac frequency does NOT follow
%                 BW-allometric scaling in neonates (HR ≈ 130–160 bpm
%                 vs. 70 bpm adult; set from clinical.HR_bpm).
%   R_shunt_pda — initial estimate from b = −1.33, but OVERRIDDEN in
%                 build_patient_params.m using Hagen-Poiseuille + clinical
%                 PDA diameter and Doppler pressure gradient.
%   V0_lv, V0_rv, V0_la, V0_ra — unstressed volumes: set from clinical
%                 body weight directly in build_patient_params.m.
%
% INPUTS:
%   params_adult  - struct of ADULT reference LPM parameters (unchanged
%                   Ortiz-Rangel 2022 values)
%   BW_neo_kg     - neonatal/patient body weight                [kg]
%   BW_adult_kg   - adult reference body weight (default 70 kg) [kg]
%
% OUTPUTS:
%   params_neo  - copy of params_adult with allometric scaling applied
%                 to all appropriate fields (see table above)
%
% USAGE EXAMPLE:
%   p_adult = default_parameters_adult();
%   p_neo   = scale_params_allometric(p_adult, 1.237, 70);
%
% REFERENCES:
%   [P&F2000] Pennati G, Fumero R. (2000). Ann Biomed Eng 28:442–452.
%   [S2026]   Seemann G et al. (2026). ASAIO J 72(3):207–215.
%   [OR2022]  Ortiz-Rangel et al. (2022). Biomed Signal Process Control
%             71:103151.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2026-05-12
% VERSION:  1.0
% -----------------------------------------------------------------------

%% Input validation
if BW_neo_kg <= 0
    error('SCALE_PARAMS_ALLOMETRIC: BW_neo_kg must be > 0. Got %.3f', BW_neo_kg);
end
if BW_adult_kg <= 0
    error('SCALE_PARAMS_ALLOMETRIC: BW_adult_kg must be > 0. Got %.3f', BW_adult_kg);
end

params_neo = params_adult;   % Deep copy; only scaled fields will change

%% -----------------------------------------------------------------------
%  Allometric ratio (dimensionless)
% -----------------------------------------------------------------------
bw_ratio = BW_neo_kg / BW_adult_kg;

fprintf('SCALE_PARAMS_ALLOMETRIC: BW ratio = %.4f (%.3f kg / %.1f kg)\n', ...
    bw_ratio, BW_neo_kg, BW_adult_kg);

%% -----------------------------------------------------------------------
%  Allometric exponents — Pennati & Fumero (2000) [P&F2000]
% -----------------------------------------------------------------------
% Viscous resistance: blood viscosity × vessel geometry → b = −1.0
%   R ∝ μ L / r^4; vessel radii and lengths scale with BW^(1/3) → b = −1
%   Reference: [P&F2000] Table 1, viscous resistance row
b_R_viscous    = -0.67;

% Turbulent/convective resistance: b = −1.33 (accounts for inertial effects)
%   Reference: [P&F2000] Table 1, turbulent resistance row
b_R_turbulent  = -1.33;

% Compliance: C ∝ vessel volume / elastic modulus → b = +1.33
%   Reference: [P&F2000] Table 1, compliance row
b_C_compliance = +1.33;

% Inertance: L ∝ ρ × length / cross-section area → b = −0.33
%   Reference: [P&F2000] Table 1, inertance row
b_L_inertance  = -0.33;

% Cardiac elastance: E scales inversely with BW (b ≈ −1.0)
%   Rationale: Ventricular pressure targets decrease (≈50–70 mmHg LV vs.
%   ≈120 mmHg adult) AND stroke volume scales ∝ BW, so E = P/SV → b ≈ −1.
%   This is consistent with Seemann et al. (2026) who report linear SV
%   scaling with weight across birth-to-adult spectrum.
%   Reference: [S2026] §2.3 cardiac scaling; [OR2022] Table 2 adult baseline
b_E_elastance  = -1.00;

%% -----------------------------------------------------------------------
%  SECTION 1: Vascular Resistances (viscous) — b = −1.0 [P&F2000]
%  Applied to: valve resistances and vascular bed resistances
%  NOT applied to: R_systemic (clinical-derived), R_shunt_pda (geometry)
% -----------------------------------------------------------------------
sf_R = bw_ratio ^ b_R_viscous;   % Scaling factor for viscous resistance

% Valve resistances
params_neo.R_tv       = params_adult.R_tv       * sf_R;  % Tricuspid valve
params_neo.R_pv_valve = params_adult.R_pv_valve * sf_R;  % Pulmonary valve
params_neo.R_mv       = params_adult.R_mv       * sf_R;  % Mitral valve
params_neo.R_av       = params_adult.R_av       * sf_R;  % Aortic valve

% Vascular bed resistances
params_neo.R_pa       = params_adult.R_pa       * sf_R;  % Pulmonary artery resistance
params_neo.R_pv_veins = params_adult.R_pv_veins * sf_R;  % Pulmonary venous resistance
params_neo.R_ao       = params_adult.R_ao       * sf_R;  % Aortic characteristic resistance

% R_systemic: NOT scaled allometrically.
% Must be set from clinical MAP/CO in build_patient_params.m.
% (See function header for rationale.)

fprintf('  [R_viscous] Scaling factor (b=%.2f): %.6f\n', b_R_viscous, sf_R);
fprintf('    R_tv       = %.5f mmHg·s/mL  (adult: %.5f)\n', params_neo.R_tv,       params_adult.R_tv);
fprintf('    R_pv_valve = %.5f mmHg·s/mL  (adult: %.5f)\n', params_neo.R_pv_valve, params_adult.R_pv_valve);
fprintf('    R_mv       = %.5f mmHg·s/mL  (adult: %.5f)\n', params_neo.R_mv,       params_adult.R_mv);
fprintf('    R_av       = %.5f mmHg·s/mL  (adult: %.5f)\n', params_neo.R_av,       params_adult.R_av);
fprintf('    R_pa       = %.5f mmHg·s/mL  (adult: %.5f)\n', params_neo.R_pa,       params_adult.R_pa);
fprintf('    R_pv_veins = %.5f mmHg·s/mL  (adult: %.5f)\n', params_neo.R_pv_veins, params_adult.R_pv_veins);
fprintf('    R_ao       = %.5f mmHg·s/mL  (adult: %.5f)\n', params_neo.R_ao,       params_adult.R_ao);

%% -----------------------------------------------------------------------
%  SECTION 2: PDA Shunt Resistance (turbulent) — b = −1.33 [P&F2000]
%  This is an INITIAL estimate using the turbulent exponent.
%  In build_patient_params.m this value is OVERRIDDEN by a Hagen-Poiseuille
%  estimate from the clinical PDA diameter and Doppler pressure gradient,
%  which guarantees Qp/Qs > 1.0 for L→R shunt.
% -----------------------------------------------------------------------
sf_R_turb = bw_ratio ^ b_R_turbulent;

params_neo.R_shunt_pda = params_adult.R_shunt_pda * sf_R_turb;
% Note: R_shunt_pda = Inf in adult baseline (no PDA); overridden below
% in build_patient_params regardless. This line handles cases where a
% non-Inf adult R_shunt_pda seed is ever provided.

%% -----------------------------------------------------------------------
%  SECTION 3: Vascular Compliances — b = +1.33 [P&F2000]
%  Applied to: all vascular bed and atrial compliances
%  C_rv and C_lv implicit via time-varying elastance (not scaled here)
% -----------------------------------------------------------------------
sf_C = bw_ratio ^ b_C_compliance;   % Scaling factor for compliance

params_neo.C_ao  = params_adult.C_ao  * sf_C;  % Aortic compliance
params_neo.C_sys = params_adult.C_sys * sf_C;  % Systemic venous compliance
params_neo.C_pa  = params_adult.C_pa  * sf_C;  % Pulmonary artery compliance
params_neo.C_pv  = params_adult.C_pv  * sf_C;  % Pulmonary venous compliance
params_neo.C_ra  = params_adult.C_ra  * sf_C;  % Right atrial compliance
params_neo.C_la  = params_adult.C_la  * sf_C;  % Left atrial compliance

fprintf('  [C_compliance] Scaling factor (b=%.2f): %.6f\n', b_C_compliance, sf_C);
fprintf('    C_ao  = %.5f mL/mmHg  (adult: %.5f)\n', params_neo.C_ao,  params_adult.C_ao);
fprintf('    C_sys = %.5f mL/mmHg  (adult: %.5f)\n', params_neo.C_sys, params_adult.C_sys);
fprintf('    C_pa  = %.5f mL/mmHg  (adult: %.5f)\n', params_neo.C_pa,  params_adult.C_pa);
fprintf('    C_pv  = %.5f mL/mmHg  (adult: %.5f)\n', params_neo.C_pv,  params_adult.C_pv);
fprintf('    C_ra  = %.5f mL/mmHg  (adult: %.5f)\n', params_neo.C_ra,  params_adult.C_ra);
fprintf('    C_la  = %.5f mL/mmHg  (adult: %.5f)\n', params_neo.C_la,  params_adult.C_la);

%% -----------------------------------------------------------------------
%  SECTION 4: Vascular Inertances — b = −0.33 [P&F2000]
%  Applied to: pulmonary and aortic inertances; PDA inertance
% -----------------------------------------------------------------------
sf_L = bw_ratio ^ b_L_inertance;   % Scaling factor for inertance

params_neo.L_pa        = params_adult.L_pa        * sf_L;  % Pulmonary artery inertance
params_neo.L_ao        = params_adult.L_ao         * sf_L;  % Aortic inertance
params_neo.L_shunt_pda = params_adult.L_shunt_pda  * sf_L;  % PDA segment inertance

fprintf('  [L_inertance] Scaling factor (b=%.2f): %.6f\n', b_L_inertance, sf_L);
fprintf('    L_pa        = %.3e mmHg·s²/mL  (adult: %.3e)\n', params_neo.L_pa,        params_adult.L_pa);
fprintf('    L_ao        = %.3e mmHg·s²/mL  (adult: %.3e)\n', params_neo.L_ao,        params_adult.L_ao);
fprintf('    L_shunt_pda = %.3e mmHg·s²/mL  (adult: %.3e)\n', params_neo.L_shunt_pda, params_adult.L_shunt_pda);

%% -----------------------------------------------------------------------
%  SECTION 5: Cardiac Elastances — b = −1.0 [S2026, OR2022]
%  Rationale:
%   Adult LV systolic pressure ≈ 120 mmHg, SV ≈ 70 mL → Emax ≈ 1.7 mmHg/mL
%   Neonate LV systolic pressure ≈ 50–60 mmHg, SV ≈ 5–7 mL → Emax ≈ 8–12 mmHg/mL
%   But: Emax_neo/Emax_adult = (P_neo/P_adult) × (SV_adult/SV_neo)
%                             ≈ (55/120) × (70/6) ≈ 0.46 × 11.7 ≈ 5.4
%   This equals (BW_neo/BW_adult)^(−1) = (1.237/70)^(-1) ≈ 56.6 ... too high.
%
%   Practical correction: The FLAT 0.8 multiplier applied in the original
%   default_parameters.m (Emax_lv = 2.0 × 0.8 = 1.6 mmHg/mL) was derived
%   for an ADULT, not a neonate. It fails because:
%     (a) Adult Emax_lv ≈ 2.0 mmHg/mL yields SV ≈ 70 mL when P_lv ≈ 120 mmHg.
%         For a neonate with P_lv ≈ 55 mmHg and SV ≈ 6 mL, the true Emax
%         should be ≈ 9 mmHg/mL — roughly 4.5× adult, not 0.8× adult.
%     (b) The ×0.8 scale was insufficient; it only mildly reduced SV.
%         The resulting SV ≈ 44 mL confirms the cardiac pump was still
%         calibrated for an adult-sized heart.
%   With b = −1.0 and BW_ratio ≈ 0.0177:
%     Emax_lv_neo = 2.0 × (0.0177)^(−1) ≈ 113 mmHg/mL (too large raw)
%   Therefore, the ABSOLUTE elastance values are SET from clinical SV and
%   MAP targets in build_patient_params.m (Emax = P_lv_sys / SV), and the
%   allometric scaling here provides an INITIAL SEED only. The clinical
%   calibration overwrites these values.
%
%   Reference: [S2026] §2.3; [OR2022] Table 2
% -----------------------------------------------------------------------
sf_E = bw_ratio ^ b_E_elastance;

params_neo.Emax_lv = params_adult.Emax_lv * sf_E;
params_neo.Emin_lv = params_adult.Emin_lv * sf_E;
params_neo.Emax_rv = params_adult.Emax_rv * sf_E;
params_neo.Emin_rv = params_adult.Emin_rv * sf_E;

% Clamp elastances to physiologically plausible neonatal range
% Emax_lv for neonate: typically 5–15 mmHg/mL (Senzaki 1996, Pettersen 2008)
% Emax_rv for neonate: typically 2–8 mmHg/mL
EMAX_LV_MIN = 3.0;   EMAX_LV_MAX = 20.0;   % [mmHg/mL]
EMAX_RV_MIN = 1.5;   EMAX_RV_MAX = 12.0;   % [mmHg/mL]

params_neo.Emax_lv = max(EMAX_LV_MIN, min(EMAX_LV_MAX, params_neo.Emax_lv));
params_neo.Emin_lv = params_neo.Emax_lv * 0.05;   % 5% of Emax (Stergiopulos 1996)
params_neo.Emax_rv = max(EMAX_RV_MIN, min(EMAX_RV_MAX, params_neo.Emax_rv));
params_neo.Emin_rv = params_neo.Emax_rv * 0.05;

fprintf('  [E_elastance] Scaling factor (b=%.2f): %.6f (clamped to physiological range)\n', ...
    b_E_elastance, sf_E);
fprintf('    Emax_lv = %.4f mmHg/mL  (adult: %.4f)\n', ...
    params_neo.Emax_lv, params_adult.Emax_lv);
fprintf('    Emin_lv = %.4f mmHg/mL  (adult: %.4f)\n', params_neo.Emin_lv, params_adult.Emin_lv);
fprintf('    Emax_rv = %.4f mmHg/mL  (adult: %.4f)\n', params_neo.Emax_rv, params_adult.Emax_rv);
fprintf('    Emin_rv = %.4f mmHg/mL  (adult: %.4f)\n', params_neo.Emin_rv, params_adult.Emin_rv);

%% -----------------------------------------------------------------------
%  SECTION 6: Unstressed volumes — scaled with BW (b = +1.0)
%  V0 ∝ BW (cardiac chamber volume scales approximately linearly with BW)
%  Reference: [S2026] §2.2
% -----------------------------------------------------------------------
sf_V = bw_ratio ^ 1.0;

params_neo.V0_lv = params_adult.V0_lv * sf_V;   % LV dead volume [mL]
params_neo.V0_rv = params_adult.V0_rv * sf_V;   % RV dead volume [mL]
params_neo.V0_la = params_adult.V0_la * sf_V;   % LA dead volume [mL]
params_neo.V0_ra = params_adult.V0_ra * sf_V;   % RA dead volume [mL]

% Clamp to minimum physiological values for a 1 kg neonate [mL]
params_neo.V0_lv = max(0.5, params_neo.V0_lv);
params_neo.V0_rv = max(0.5, params_neo.V0_rv);
params_neo.V0_la = max(0.3, params_neo.V0_la);
params_neo.V0_ra = max(0.3, params_neo.V0_ra);

fprintf('  [V0_volumes] Scaling factor (b=1.00): %.6f\n', sf_V);
fprintf('    V0_lv = %.4f mL  |  V0_rv = %.4f mL\n', params_neo.V0_lv, params_neo.V0_rv);
fprintf('    V0_la = %.4f mL  |  V0_ra = %.4f mL\n', params_neo.V0_la, params_neo.V0_ra);

%% -----------------------------------------------------------------------
%  Parameters explicitly NOT scaled (document exclusions)
% -----------------------------------------------------------------------
% params_neo.R_systemic     ← NOT scaled; derived from clinical MAP/CO
% params_neo.HR_bpm         ← NOT scaled; set from clinical.HR_bpm
% params_neo.T_cardiac      ← NOT scaled; derived from HR
% params_neo.Ts1, .Ts2      ← NOT scaled; derived from T_cardiac
% params_neo.R_shunt_pda    ← overridden in build_patient_params (geometry)
% params_neo.rho_blood_kg_m3 ← physical constant; unchanged
% params_neo.mu_blood_Pa_s   ← physical constant; unchanged
% CoA params (R_coa_viscous, K_turb_coa, L_coa) ← geometry-derived in
%   build_coa_params.m; not allometrically scaled

fprintf('  [EXCLUDED from allometric scaling]: R_systemic, HR, R_shunt_pda (overridden),\n');
fprintf('    rho_blood, mu_blood, CoA params (geometry-derived)\n\n');

end
