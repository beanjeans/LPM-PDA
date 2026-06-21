function params_coa = build_coa_params(params_pda, clinical, stenosis_pct, coa_length_mm)
% BUILD_COA_PARAMS
% -----------------------------------------------------------------------
% Extends the PDA parameter struct with a virtual Coarctation of the Aorta
% (CoA) at a specified stenosis severity and anatomical length.
%
% CoA geometry is derived from the patient's ascending aorta diameter
% (D_aao_mm) as the reference normal diameter, then stenosed:
%
%   D_coa_eff = D_ref * sqrt(1 - stenosis_pct/100)
%
% INPUTS:
%   params_pda    - PDA-calibrated parameter struct (from build_patient_params)
%   clinical      - clinical struct (from load_patient_data)
%   stenosis_pct  - virtual CoA anatomical severity         [%], range (0,100)
%   coa_length_mm - CoA segment length                      [mm]
%                   If omitted or empty, defaults to 3 mm (discrete CoA).
%                   Typical scenario values:
%                     3 mm  — discrete/short CoA  (< 5 mm threshold)
%                     5 mm  — threshold case       (boundary)
%                     8 mm  — long-segment CoA
%                    10 mm  — severe long-segment CoA
%
% OUTPUTS:
%   params_coa    - extended parameter struct with CoA fields added and
%                   scenario_coa = true.  Key new fields:
%       .coa_length_mm          — actual CoA length used           [mm]
%       .coa_length_category    — 'discrete/short' or 'long-segment'
%       .severity_thresholds    — editable struct with DP thresholds [mmHg]
%       .R_coa_viscous          — viscous Poiseuille resistance [mmHg·s/mL]
%       .K_turb_coa             — Bernoulli turbulent coeff    [mmHg/(mL/s)²]
%       .L_coa                  — segment inertance             [mmHg·s²/mL]
%       .D_coa_eff_mm           — effective stenotic diameter        [mm]
%       .D_ref_mm               — reference (normal) aortic diameter [mm]
%
% SEVERITY THRESHOLDS (edit params_coa.severity_thresholds to change):
%   The predicted_CoA_severity reported by compute_clinical_indices is
%   based on the simulated pressure gradient — NOT on stenosis_pct directly.
%   Thresholds follow ESC/ACC guideline consensus (Ref [3]):
%     Mild:       ΔP_mean_systolic  <  20 mmHg
%     Moderate:   ΔP_mean_systolic  20–40 mmHg
%     Severe:     ΔP_mean_systolic  >  40 mmHg
%
% CoA LENGTH CATEGORIES (Ref [4]):
%   Discrete/short CoA:  length < 5 mm
%   Long-segment CoA:    length ≥ 5 mm
%
% PDA ROLE:
%   PDA parameters (R_shunt_pda, L_shunt_pda) are RETAINED in params_coa.
%   A patent PDA may raise distal aortic pressure (P_ao_dist) via retrograde
%   flow, thereby reducing the observed ΔP_CoA and MASKING true CoA severity.
%   PDA size must therefore be reported as a haemodynamic MODIFIER alongside
%   the ΔP_CoA-based severity classification.
%
% SIGN CONVENTIONS:
%   - Q_coa > 0        : forward flow, proximal → distal aorta
%   - DeltaP_coa       = P_ao (proximal) − P_ao_dist (distal)
%   - Q_shunt_pda > 0  : left-to-right PDA shunt (Ao → PA)
%
% REFERENCES:
%   [1] Keshavarz-Motamed et al. (2011). J Biomech 44:2817–2825.
%       (Viscous + turbulent CoA model, Eq. 3–5)
%   [2] Vergales et al. (2013). Pediatr Cardiol 34:1616–1623.
%       (Native neonatal coarctation anatomy)
%   [3] Baumgartner et al. (2010). Eur Heart J 31(19):2369–2417.
%       (ESC valvular/CoA guidelines: significant gradient ≥ 20 mmHg;
%        severe gradient > 40 mmHg — peak-to-peak/Doppler)
%   [4] Vergales et al. (2013) / Campbell et al. (2002).
%       (Length classification: discrete < 5 mm; long-segment ≥ 5 mm)
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  2.0  — pressure-gradient severity classification; variable length
% -----------------------------------------------------------------------

uc          = unit_conversion();
params_coa  = params_pda;     % Start from PDA-calibrated params

%% -----------------------------------------------------------------------
%  0. Handle optional coa_length_mm argument
% -----------------------------------------------------------------------
if nargin < 4 || isempty(coa_length_mm)
    coa_length_mm = 3;   % Default: discrete CoA scenario [mm]
    fprintf('  [build_coa_params] coa_length_mm not supplied → using %.0f mm (discrete CoA default)\n', ...
            coa_length_mm);
end

%% -----------------------------------------------------------------------
%  1. Validate inputs
% -----------------------------------------------------------------------
if stenosis_pct <= 0 || stenosis_pct >= 100
    error('BUILD_COA_PARAMS: stenosis_pct must be in (0, 100). Got %.1f', stenosis_pct);
end
if coa_length_mm <= 0
    error('BUILD_COA_PARAMS: coa_length_mm must be > 0. Got %.2f', coa_length_mm);
end

%% -----------------------------------------------------------------------
%  2. CoA Length Category
%  Literature threshold: discrete < 5 mm; long-segment ≥ 5 mm — Ref [4]
% -----------------------------------------------------------------------
COA_LENGTH_THRESHOLD_MM = 5;   % [mm] — edit here to change category boundary

params_coa.coa_length_mm = coa_length_mm;   % [mm] — store as-supplied

if coa_length_mm < COA_LENGTH_THRESHOLD_MM
    params_coa.coa_length_category = 'discrete/short';
else
    params_coa.coa_length_category = 'long-segment';
end

L_coa_m = coa_length_mm * uc.mm_to_m;      % Convert to [m] for physics

%% -----------------------------------------------------------------------
%  3. Pressure Gradient Severity Thresholds
%  These govern predicted_CoA_severity in compute_clinical_indices.
%  Edit this struct only — nowhere else in the codebase — Ref [3]
% -----------------------------------------------------------------------
sev_thr.mild_upper_mmHg     = 20;   % ΔP < 20 mmHg → Mild
sev_thr.moderate_upper_mmHg = 40;   % 20–40 mmHg → Moderate
%                                    % ΔP > 40 mmHg → Severe/Critical

params_coa.severity_thresholds = sev_thr;

%% -----------------------------------------------------------------------
%  4. CoA Geometry (from stenosis percentage)
% -----------------------------------------------------------------------
% Reference aorta diameter = ascending aorta from echocardiography
D_ref_m     = clinical.D_aao_mm * uc.mm_to_m;       % [mm] → [m]

% Effective CoA orifice diameter
%   Area stenosis: A_coa / A_ref = 1 - stenosis_pct/100
%   Diameter:      D_coa = D_ref * sqrt(1 - stenosis_pct/100)
D_coa_eff_m = D_ref_m * sqrt(1 - stenosis_pct / 100);  % [m]
A_coa_m2    = pi * (D_coa_eff_m / 2)^2;                 % [m²]

%% -----------------------------------------------------------------------
%  5. Viscous Resistance (Poiseuille flow in stenotic segment)
%  R_viscous = (128 × μ × L) / (π × D⁴)    [Pa·s/m³]  — Ref [1] Eq. 3
% -----------------------------------------------------------------------
R_coa_viscous_SI = (128 * params_coa.mu_blood_Pa_s * L_coa_m) ...
                 / (pi * D_coa_eff_m^4);               % [Pa·s/m³]

% Convert to model units [mmHg·s/mL]
params_coa.R_coa_viscous = R_coa_viscous_SI * uc.Pa_s_m3_to_mmHg_s_mL;

%% -----------------------------------------------------------------------
%  6. Turbulent (Bernoulli) Coefficient — Ref [1] Eq. 4–5
%  ΔP_turb = 0.5 × ρ × (Q/A)²  →  ΔP = K_turb × Q × |Q|
%  K_turb_SI = 0.5 × ρ / A²    [Pa / (m³/s)²]
% -----------------------------------------------------------------------
K_turb_SI            = 0.5 * params_coa.rho_blood_kg_m3 / (A_coa_m2^2);  % [Pa/(m³/s)²]
params_coa.K_turb_coa = K_turb_SI * uc.Pa_m3s2_to_mmHg_mLs2;             % [mmHg/(mL/s)²]

%% -----------------------------------------------------------------------
%  7. CoA Segment Inertance — Ref [1]
%  L_inertance = (4 × ρ × L) / (π × D²)   [Pa·s²/m³]
% -----------------------------------------------------------------------
L_coa_SI         = (4 * params_coa.rho_blood_kg_m3 * L_coa_m) ...
                 / (pi * D_coa_eff_m^2);               % [Pa·s²/m³]
params_coa.L_coa = L_coa_SI * uc.Pa_s2_m3_to_mmHg_s2_mL;                 % [mmHg·s²/mL]

%% -----------------------------------------------------------------------
%  8. Extend State Vector for CoA
%  New states: P_ao_dist (distal aortic pressure) and Q_coa (CoA flow)
% -----------------------------------------------------------------------
params_coa.idx.P_ao_dist = 12;  % Distal aortic pressure     [mmHg]
params_coa.idx.Q_coa     = 13;  % Flow through CoA stenosis  [mL/s]

% Initial conditions: distal Ao starts ~15% below proximal; CoA flow = 0
params_coa.X0(12) = params_coa.X0(params_coa.idx.P_ao) * 0.85;  % [mmHg]
params_coa.X0(13) = 0;                                            % [mL/s]

%% -----------------------------------------------------------------------
%  9. Set Scenario Flag and Metadata
% -----------------------------------------------------------------------
params_coa.scenario_coa   = true;
params_coa.stenosis_pct   = stenosis_pct;                        % [%]
params_coa.D_coa_eff_mm   = D_coa_eff_m * uc.m_to_mm;           % [mm] — for reporting
params_coa.D_ref_mm       = clinical.D_aao_mm;                   % [mm]

%% -----------------------------------------------------------------------
%  10. Display CoA parameters
% -----------------------------------------------------------------------
fprintf('--- VIRTUAL CoA PARAMETERS (Stenosis: %.0f%% | Length: %.1f mm — %s) ---\n', ...
    stenosis_pct, coa_length_mm, params_coa.coa_length_category);
fprintf('  D_ref (AscAo):     %.2f mm\n', params_coa.D_ref_mm);
fprintf('  D_coa_eff:         %.2f mm\n', params_coa.D_coa_eff_mm);
fprintf('  CoA length:        %.1f mm  [%s CoA]\n', coa_length_mm, params_coa.coa_length_category);
fprintf('  R_coa_viscous:     %.4f mmHg·s/mL\n', params_coa.R_coa_viscous);
fprintf('  K_turb_coa:        %.6f mmHg/(mL/s)²\n', params_coa.K_turb_coa);
fprintf('  L_coa:             %.2e mmHg·s²/mL\n', params_coa.L_coa);
fprintf('  Severity thresholds: Mild <%.0f mmHg | Moderate %.0f–%.0f mmHg | Severe >%.0f mmHg\n\n', ...
    sev_thr.mild_upper_mmHg, sev_thr.mild_upper_mmHg, ...
    sev_thr.moderate_upper_mmHg, sev_thr.moderate_upper_mmHg);

end
