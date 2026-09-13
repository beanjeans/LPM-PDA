function params = default_parameters(BW_neo_kg)
% DEFAULT_PARAMETERS
% -----------------------------------------------------------------------
% Returns the reference parameter set for a neonatal cardiovascular LPM,
% scaled allometrically from adult reference values using Pennati & Fumero
% (2000) exponents applied per parameter type.
%
% SCALING ARCHITECTURE:
%   1. Adult reference values (Ortiz-Rangel 2022, Table 2) are stored
%      UNCHANGED in params.adult_ref.* — never hardcode them elsewhere.
%   2. scale_params_allometric() converts adult → neonatal using:
%      - Compliance:          b = +1.33  [Pennati & Fumero 2000]
%      - Viscous resistance:  b = −1.00  [Pennati & Fumero 2000]
%      - Inertance:           b = −0.33  [Pennati & Fumero 2000]
%      - Cardiac elastance:   b = −1.00  [Seemann et al. 2026]
%   3. R_systemic, HR and PDA resistance are NOT allometrically scaled —
%      they are derived from clinical data in build_patient_params.m.
%
% INPUTS (required):
%   BW_neo_kg  - patient body weight [kg], converted from CSV column BB (grams).
%                e.g.: BW_neo_kg = clinical.weight_g / 1000
%                NO default — must always be supplied from clinical data.
%
% OUTPUTS:
%   params       - struct of all LPM parameters for neonatal patient,
%                  including state index map, initial conditions,
%                  and params.adult_ref with unmodified adult values.
%
% REFERENCES:
%   [OR2022]  Ortiz-Rangel et al. (2022). Biomed Signal Process Control
%             71:103151.  (Adult baseline Table 2)
%   [P&F2000] Pennati G, Fumero R. (2000). Ann Biomed Eng 28:442–452.
%             DOI: 10.1114/1.282.  (Allometric exponents per param type)
%   [S2026]   Seemann G et al. (2026). ASAIO J 72(3):207–215.
%             DOI: 10.1097/MAT.0000000000002528.  (Birth-to-adult LPM
%             allometric validation; elastance scaling rationale)
%   [Rud2001] Rudolph AM (2001). Congenital Diseases of the Heart.
%             (Neonatal haemodynamic reference values)
%   [Ste1996] Stergiopulos et al. (1996). Am J Physiol 270(6):H2050–H2059.
%             (Elastance normalisation; Ts1/Ts2 timing)
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2026-05-12
% VERSION:  2.0  — allometric scaling replaces flat multipliers
% -----------------------------------------------------------------------

BW_ADULT_KG = 70;        % [kg] — Ortiz-Rangel adult reference body weight [OR2022]

%% -----------------------------------------------------------------------
%  STATE VECTOR INDEX STRUCT
%  Units: pressures [mmHg], flows [mL/s]
%  Do NOT hardcode indices elsewhere (Guardrail §7.1)
% -----------------------------------------------------------------------
idx.P_ra        =  1;   % Right atrial pressure               [mmHg]
idx.P_rv        =  2;   % Right ventricular pressure          [mmHg]
idx.P_pa        =  3;   % Pulmonary artery pressure           [mmHg]
idx.Q_pa_pul    =  4;   % Pulmonary arterial flow             [mL/s]
idx.P_pv        =  5;   % Pulmonary venous pressure           [mmHg]
idx.P_la        =  6;   % Left atrial pressure                [mmHg]
idx.P_lv        =  7;   % Left ventricular pressure           [mmHg]
idx.P_ao        =  8;   % Aortic pressure (proximal)          [mmHg]
idx.Q_ao_sys    =  9;   % Systemic aortic flow                [mL/s]
idx.P_sys       = 10;   % Systemic venous pressure            [mmHg]
idx.Q_shunt_pda = 11;   % PDA shunt flow (L→R positive)      [mL/s]

%% -----------------------------------------------------------------------
%  ADULT REFERENCE PARAMETERS — Ortiz-Rangel et al. (2022) Table 2
%  These values are NEVER modified. Allometric scaling is applied in the
%  next section. Store here for traceability and reproducibility.
%
%  IMPORTANT: If you update these values, update the citation too.
% -----------------------------------------------------------------------
ar = struct();  % adult_ref namespace

% Physical constants — unchanged across body sizes
ar.rho_blood_kg_m3  = 1060;     % Blood density          [kg/m³]  — [Rud2001]
ar.mu_blood_Pa_s    = 0.004;    % Dynamic viscosity      [Pa·s]   — [Rud2001]

% Cardiac timing at adult reference HR
ar.HR_bpm           = 75;       % Adult reference HR     [bpm]    — [OR2022]
ar.T_cardiac        = 60 / ar.HR_bpm;
ar.Ts1              = 0.3  * sqrt(ar.T_cardiac);  % [s] — [Ste1996]
ar.Ts2              = 0.45 * sqrt(ar.T_cardiac);  % [s] — [Ste1996]

% Ventricular elastance — [OR2022] Table 2
% Adult Emax_lv = 2.0 mmHg/mL → gives SV ≈ 70 mL at P_lv_sys ≈ 120 mmHg
ar.Emax_lv          = 2.0;      % LV peak (end-systolic) elastance  [mmHg/mL]
ar.Emin_lv          = 0.06;     % LV diastolic elastance            [mmHg/mL]
ar.Emax_rv          = 1.0;      % RV peak elastance                 [mmHg/mL]
ar.Emin_rv          = 0.04;     % RV diastolic elastance            [mmHg/mL]

% Unstressed volumes — [OR2022] Table 2 (adult scale)
ar.V0_lv            = 5;        % LV dead volume                    [mL]
ar.V0_rv            = 5;        % RV dead volume                    [mL]
ar.V0_la            = 3;        % LA dead volume                    [mL]
ar.V0_ra            = 3;        % RA dead volume                    [mL]

% Vascular resistances [mmHg·s/mL] — [OR2022] Table 2
ar.R_tv             = 0.0025;    % Tricuspid valve resistance
ar.R_pv_valve       = 0.001;    % Pulmonary valve resistance
ar.R_mv             = 0.0025;    % Mitral valve resistance
ar.R_av             = 0.001;    % Aortic valve resistance
ar.R_pa             = 0.05;     % Pulmonary artery resistance
ar.R_pv_veins       = 0.05;     % Pulmonary venous resistance
ar.R_ao             = 0.05;     % Aortic characteristic resistance
% ar.R_systemic: NOT stored here; derived from clinical MAP/CO per patient

% Vascular compliances [mL/mmHg] — [OR2022] Table 2 adult values
%   These are the UNSCALED adult values used as the allometric base.
%   Old code used neonate_C_scale = 0.15 (flat multiplier — now replaced).
ar.C_ao             = 0.08;     % Aortic compliance
ar.C_sys            = 1.4465;     % Systemic venous compliance
ar.C_pa             = 2.17;     % Pulmonary artery compliance
ar.C_pv             = 5.0;      % Pulmonary venous compliance
ar.C_ra             = 62.42;    % Right atrial compliance
ar.C_la             = 47.65;    % Left atrial compliance

% Vascular inertances [mmHg·s²/mL] — [OR2022] Table 2
ar.L_pa             = 5e-4;     % Pulmonary artery inertance
ar.L_ao             = 5e-4;     % Aortic inertance
ar.L_shunt_pda      = 5e-4;     % PDA inertance (adult placeholder)

% PDA: absent in adult baseline
ar.R_shunt_pda      = Inf;      % No PDA in adult reference

% CoA: absent in baseline
ar.R_coa_viscous    = 0;
ar.K_turb_coa       = 0;
ar.L_coa            = 5e-5;
ar.stenosis_pct     = 0;
ar.scenario_coa     = false;

%% -----------------------------------------------------------------------
%  ALLOMETRIC SCALING: Adult → Neonatal
%  Apply Pennati & Fumero (2000) exponents via scale_params_allometric()
%  Exponent summary (see scale_params_allometric.m for full rationale):
%    Viscous resistance: b = −1.00  [P&F2000]
%    Compliance:         b = +1.33  [P&F2000]
%    Inertance:          b = −0.33  [P&F2000]
%    Cardiac elastance:  b = −1.00  [S2026]
%    Unstressed volumes: b = +1.00  [S2026]
% -----------------------------------------------------------------------
ar.idx = idx;  % Pass index map through to scaler

params_neo_scaled = scale_params_allometric(ar, BW_neo_kg, BW_ADULT_KG);

%% -----------------------------------------------------------------------
%  BUILD FINAL PARAMS STRUCT
%  Copy allometrically-scaled neonatal values; add metadata fields
% -----------------------------------------------------------------------
params = params_neo_scaled;   % Start from scaled values

% Preserve adult reference and scaling metadata
params.adult_ref        = ar;
params.BW_neo_kg        = BW_neo_kg;
params.BW_adult_kg      = BW_ADULT_KG;
params.scaling_method   = 'allometric_pennati_fumero_2000';

% Physical constants (unchanged — copy explicitly for ODE access)
params.rho_blood_kg_m3  = ar.rho_blood_kg_m3;
params.mu_blood_Pa_s    = ar.mu_blood_Pa_s;

% Index map
params.idx = idx;

%% -----------------------------------------------------------------------
%  CARDIAC TIMING — Neonatal reference HR
%  Will be overridden per patient in build_patient_params.m.
%  HR does NOT follow allometric scaling for neonates (HR ≈ 130–158 bpm
%  vs. 70 bpm adult; physiologically independent of BW in this range).
% -----------------------------------------------------------------------
HR_ref           = 130;                        % Neonatal reference HR [bpm]
params.HR_bpm    = HR_ref;
params.T_cardiac = 60 / HR_ref;               % Cycle period [s]
params.Ts1       = 0.3  * sqrt(params.T_cardiac);  % [s] — [Ste1996]
params.Ts2       = 0.45 * sqrt(params.T_cardiac);  % [s] — [Ste1996]

%% -----------------------------------------------------------------------
%  PDA PARAMETERS — default = absent; will be set in build_patient_params
% -----------------------------------------------------------------------
params.R_shunt_pda  = Inf;       % No shunt (overridden when PDA present)
% L_shunt_pda already set by allometric scaling above

%% -----------------------------------------------------------------------
%  CoA PARAMETERS — default = absent; added in build_coa_params
% -----------------------------------------------------------------------
params.scenario_coa         = false;
params.R_coa_viscous        = 0;
params.K_turb_coa           = 0;
params.L_coa                = 5e-5;     % [mmHg·s²/mL]
params.stenosis_pct         = 0;

%% -----------------------------------------------------------------------
%  INITIAL CONDITIONS [mmHg for pressures, mL/s for flows]
%  Physiologically motivated for a preterm neonate at rest — [Rud2001]
%  These are warm-start values; solver will converge to limit cycle.
% -----------------------------------------------------------------------
n_states                    = 11;
X0                          = zeros(n_states, 1);

X0(idx.P_ra)        =  4;   % [mmHg] — RA mean pressure, neonate
X0(idx.P_rv)        =  5;   % [mmHg] — RV diastolic, neonate
X0(idx.P_pa)        = 15;   % [mmHg] — PA mean, neonate (elevated PVR at birth)
X0(idx.Q_pa_pul)    =  0;   % [mL/s] — initial pulmonary flow
X0(idx.P_pv)        =  8;   % [mmHg] — pulmonary venous pressure
X0(idx.P_la)        =  5;   % [mmHg] — LA mean pressure
X0(idx.P_lv)        =  5;   % [mmHg] — LV diastolic pressure
X0(idx.P_ao)        = 42;   % [mmHg] — neonatal MAP estimate (replaces 50 mmHg adult seed)
X0(idx.Q_ao_sys)    =  0;   % [mL/s] — initial systemic flow
X0(idx.P_sys)       = 42;   % [mmHg] — systemic venous pressure seed
X0(idx.Q_shunt_pda) =  0;   % [mL/s] — initial PDA shunt flow

params.X0 = X0;

end
