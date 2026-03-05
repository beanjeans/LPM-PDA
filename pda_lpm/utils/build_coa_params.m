function params_coa = build_coa_params(params_pda, clinical, stenosis_pct)
% BUILD_COA_PARAMS
% -----------------------------------------------------------------------
% Extends the PDA parameter struct with a virtual Coarctation of the Aorta
% (CoA) at a specified stenosis severity. The CoA is modelled as a combined
% viscous (Poiseuille) and turbulent (Bernoulli) resistance in the aortic
% segment, plus an inertance term for pulsatile momentum.
%
% CoA geometry is derived from the patient's ascending aorta diameter
% (D_aao_mm) as the reference normal diameter, then stenosed:
%
%   D_coa_eff = D_ref * sqrt(1 - stenosis_pct/100)
%
% INPUTS:
%   params_pda    - PDA-calibrated parameter struct (from build_patient_params)
%   clinical      - clinical struct (from load_patient_data)
%   stenosis_pct  - virtual CoA severity          [%], in range (0, 100)
%
% OUTPUTS:
%   params_coa    - extended parameter struct with CoA fields added
%                   and scenario_coa = true
%
% ASSUMPTIONS:
%   - CoA segment length L_coa = 8 mm (typical native coarctation — Ref [2])
%   - Reference aortic diameter taken as D_aao (ascending aorta from echo)
%   - Blood: ρ = 1060 kg/m³, μ = 0.004 Pa·s
%   - State vector extended: idx.P_ao_dist (12), idx.Q_coa (13)
%
% SIGN CONVENTIONS:
%   - Q_coa > 0 : flow from proximal to distal aorta (forward)
%   - DeltaP_coa = P_ao - P_ao_dist (upstream minus downstream)
%
% REFERENCES:
%   [1] Keshavarz-Motamed et al. (2011). J Biomech 44:2817–2825.
%       (Viscous + turbulent CoA model, Eq. 3–5)
%   [2] Vergales et al. (2013). Pediatr Cardiol 34:1616–1623.
%       (L_coa ≈ 8 mm for native neonatal coarctation)
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

uc          = unit_conversion();
params_coa  = params_pda;     % Start from PDA-calibrated params

%% Validate stenosis input
if stenosis_pct <= 0 || stenosis_pct >= 100
    error('BUILD_COA_PARAMS: stenosis_pct must be in (0, 100). Got %.1f', stenosis_pct);
end

%% 1. CoA Geometry
L_coa_m       = 8e-3;                              % CoA segment length [m]  — Ref [2]

% Reference aorta diameter = ascending aorta (most reliable echo measurement)
D_ref_m       = clinical.D_aao_mm * uc.mm_to_m;  % [mm] → [m]

% Effective CoA orifice diameter at the given stenosis severity
% Area stenosis: A_coa/A_ref = 1 - stenosis_pct/100
% Diameter: D_coa = D_ref * sqrt(1 - stenosis_pct/100)
D_coa_eff_m   = D_ref_m * sqrt(1 - stenosis_pct / 100);  % [m]
A_coa_m2      = pi * (D_coa_eff_m / 2)^2;                 % [m²]

%% 2. Viscous Resistance (Poiseuille flow in stenotic segment)
% R_viscous = (128 * μ * L) / (π * D⁴)    [Pa·s/m³]
% Reference: Ref [1] Eq. 3
R_coa_viscous_SI  = (128 * params_coa.mu_blood_Pa_s * L_coa_m) ...
                  / (pi * D_coa_eff_m^4);                 % [Pa·s/m³]

% Convert to model units [mmHg·s/mL]
params_coa.R_coa_viscous = R_coa_viscous_SI * uc.Pa_s_m3_to_mmHg_s_mL;  % [mmHg·s/mL]

%% 3. Turbulent (Bernoulli) Coefficient
% DeltaP_turb = 0.5 * ρ * (Q/A)²  →  DeltaP = K_turb * Q * |Q|
% K_turb_SI = 0.5 * ρ / A²    [Pa / (m³/s)²]
% Reference: Ref [1] Eq. 4–5
K_turb_SI             = 0.5 * params_coa.rho_blood_kg_m3 / (A_coa_m2^2);  % [Pa/(m³/s)²]
params_coa.K_turb_coa = K_turb_SI * uc.Pa_m3s2_to_mmHg_mLs2;              % [mmHg/(mL/s)²]

%% 4. CoA Segment Inertance
% L_inertance = (4 * ρ * L) / (π * D²)   [Pa·s²/m³]  — Ref [1]
L_coa_SI            = (4 * params_coa.rho_blood_kg_m3 * L_coa_m) ...
                    / (pi * D_coa_eff_m^2);              % [Pa·s²/m³]
params_coa.L_coa    = L_coa_SI * uc.Pa_s2_m3_to_mmHg_s2_mL;               % [mmHg·s²/mL]

%% 5. Extend State Vector for CoA
% Add two new states: P_ao_dist (distal aortic pressure) and Q_coa (CoA flow)
params_coa.idx.P_ao_dist = 12;  % Distal aortic pressure     [mmHg]
params_coa.idx.Q_coa     = 13;  % Flow through CoA stenosis  [mL/s]

% Extend X0 for new states
params_coa.X0(12) = params_coa.X0(params_coa.idx.P_ao) * 0.85;  % [mmHg] — ~15% drop
params_coa.X0(13) = 0;                                            % [mL/s]

%% 6. Set Scenario Flag and Metadata
params_coa.scenario_coa   = true;
params_coa.stenosis_pct   = stenosis_pct;           % [%]
params_coa.D_coa_eff_mm   = D_coa_eff_m * uc.m_to_mm; % [mm] — for reporting
params_coa.D_ref_mm       = clinical.D_aao_mm;       % [mm]

%% 7. Display CoA parameters
fprintf('--- VIRTUAL CoA PARAMETERS (Stenosis: %.0f%%) ---\n', stenosis_pct);
fprintf('  D_ref (AscAo):   %.2f mm\n', params_coa.D_ref_mm);
fprintf('  D_coa_eff:       %.2f mm\n', params_coa.D_coa_eff_mm);
fprintf('  R_coa_viscous:   %.4f mmHg·s/mL\n', params_coa.R_coa_viscous);
fprintf('  K_turb_coa:      %.6f mmHg/(mL/s)²\n', params_coa.K_turb_coa);
fprintf('  L_coa:           %.2e mmHg·s²/mL\n\n', params_coa.L_coa);

end
