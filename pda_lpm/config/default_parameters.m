function params = default_parameters()
% DEFAULT_PARAMETERS
% -----------------------------------------------------------------------
% Returns the reference parameter set for a healthy neonate cardiovascular
% LPM. All values are sourced from published literature. Patient-specific
% overrides are applied in build_patient_params.m.
%
% OUTPUTS:
%   params  - struct of all LPM parameters (see sections below)
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%       (Table 2 — healthy adult; scaled here for neonate)
%   [2] Stergiopulos et al. (1996). Am J Physiol 270(6):H2050–H2059.
%       (Elastance normalisation)
%   [3] Rudolph AM (2001). Congenital Diseases of the Heart. Futura Pub.
%       (Neonate haemodynamic reference values)
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

%% -----------------------------------------------------------------------
%  STATE VECTOR INDEX STRUCT
%  Define once here; never hardcode indices elsewhere (Guardrail §7.1)
%  Units: pressures [mmHg], flows [mL/s]
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
% PDA extension (added dynamically in build_patient_params)
idx.Q_shunt_pda = 11;   % PDA shunt flow (L→R positive)      [mL/s]

params.idx = idx;

%% -----------------------------------------------------------------------
%  PHYSICAL CONSTANTS
% -----------------------------------------------------------------------
params.rho_blood_kg_m3  = 1060;   % Blood density          [kg/m³]  — Ref [3]
params.mu_blood_Pa_s    = 0.004;  % Dynamic viscosity      [Pa·s]   — Ref [3]

%% -----------------------------------------------------------------------
%  CARDIAC TIMING (scaled for neonate HR ≈ 130 bpm)
%  Will be overridden per patient in build_patient_params.m
% -----------------------------------------------------------------------
HR_ref              = 130;                % Reference HR        [bpm]
params.HR_bpm       = HR_ref;
params.T_cardiac    = 60 / HR_ref;        % Cycle period        [s]
params.Ts1          = 0.3 * sqrt(params.T_cardiac);  % Contraction onset [s] — Ref [2]
params.Ts2          = 0.45 * sqrt(params.T_cardiac); % Relaxation end    [s] — Ref [2]

%% -----------------------------------------------------------------------
%  VENTRICULAR ELASTANCE (Neonate scaled from Ref [1] Table 2)
%  Adult: Emax_lv = 2.0, Emin_lv = 0.06 — scaled ×0.8 for smaller SV
% -----------------------------------------------------------------------
params.Emax_lv      = 1.6;     % LV peak (end-systolic) elastance  [mmHg/mL] — Ref [1] ×0.8
params.Emin_lv      = 0.06;    % LV diastolic elastance            [mmHg/mL] — Ref [1]
params.Emax_rv      = 0.8;     % RV peak elastance                 [mmHg/mL] — Ref [1] ×0.8
params.Emin_rv      = 0.04;    % RV diastolic elastance            [mmHg/mL] — Ref [1]

% Unstressed (zero-pressure) volumes — Refs [1,3]
params.V0_lv        = 5;       % LV dead volume                    [mL]
params.V0_rv        = 5;       % RV dead volume                    [mL]
params.V0_la        = 3;       % LA dead volume                    [mL]
params.V0_ra        = 3;       % RA dead volume                    [mL]

%% -----------------------------------------------------------------------
%  VASCULAR RESISTANCES  [mmHg·s/mL]
%  Adult values from Ref [1] Table 2; neonate scaling factor ≈ 8–10×
%  (smaller vessels → higher resistance per unit body size)
% -----------------------------------------------------------------------
% Valve resistances (neonate thin leaflets — low resistance)
params.R_tv         = 0.005;   % Tricuspid valve                   [mmHg·s/mL] — Ref [1]
params.R_pv_valve   = 0.005;   % Pulmonary valve                   [mmHg·s/mL] — Ref [1]
params.R_mv         = 0.005;   % Mitral valve                      [mmHg·s/mL] — Ref [1]
params.R_av         = 0.005;   % Aortic valve                      [mmHg·s/mL] — Ref [1]

% Vascular bed resistances
params.R_pa         = 0.05;    % Pulmonary artery resistance        [mmHg·s/mL] — Ref [1]
params.R_pv_veins   = 0.05;    % Pulmonary venous resistance        [mmHg·s/mL] — Ref [1]
params.R_ao         = 0.05;    % Aortic characteristic resistance   [mmHg·s/mL] — Ref [1]
% R_systemic: computed from MAP/CO in build_patient_params — not hardcoded

%% -----------------------------------------------------------------------
%  VASCULAR COMPLIANCES  [mL/mmHg]
%  Neonate scaling factor ≈ 0.15× adult (body weight ~3 kg vs ~70 kg)
%  Adult values from Ref [1] Table 2; neonate: multiply by 0.15
% -----------------------------------------------------------------------
neonate_C_scale     = 0.15;    % Neonate compliance scaling factor [dimensionless] — Ref [3]

params.C_ao         = 0.08  * neonate_C_scale;  % Aortic compliance          [mL/mmHg]
params.C_sys        = 1.45  * neonate_C_scale;  % Systemic venous compliance [mL/mmHg]
params.C_pa         = 2.17  * neonate_C_scale;  % Pulmonary artery compliance[mL/mmHg]
params.C_pv         = 5.0   * neonate_C_scale;  % Pulmonary venous compliance[mL/mmHg]
params.C_ra         = 62.42 * neonate_C_scale;  % Right atrial compliance    [mL/mmHg]
params.C_la         = 47.65 * neonate_C_scale;  % Left atrial compliance     [mL/mmHg]
% C_rv and C_lv are implicitly represented via time-varying elastance

%% -----------------------------------------------------------------------
%  VASCULAR INERTANCES  [mmHg·s²/mL]
%  From Ref [1] Table 2 — inertance scales weakly with body size
% -----------------------------------------------------------------------
params.L_pa         = 5e-5;    % Pulmonary artery inertance        [mmHg·s²/mL] — Ref [1]
params.L_ao         = 5e-5;    % Aortic inertance                  [mmHg·s²/mL] — Ref [1]

%% -----------------------------------------------------------------------
%  PDA PARAMETERS (default = absent; set in build_patient_params)
%  R_shunt_pda = Inf means no shunt
% -----------------------------------------------------------------------
params.R_shunt_pda  = Inf;     % PDA resistance                    [mmHg·s/mL]
params.L_shunt_pda  = 5e-5;    % PDA inertance                     [mmHg·s²/mL] — Ref [1]

%% -----------------------------------------------------------------------
%  CoA PARAMETERS (default = absent; added in build_coa_params)
% -----------------------------------------------------------------------
params.scenario_coa         = false;    % CoA scenario flag
params.R_coa_viscous        = 0;        % Viscous CoA resistance    [mmHg·s/mL]
params.K_turb_coa           = 0;        % Turbulent coefficient     [mmHg/(mL/s)²]
params.L_coa                = 5e-5;     % CoA segment inertance     [mmHg·s²/mL]
params.stenosis_pct         = 0;        % Stenosis severity         [%]

%% -----------------------------------------------------------------------
%  INITIAL CONDITIONS  [mmHg for pressures, mL/s for flows]
%  Physiologically motivated for a neonate at rest — Ref [3]
% -----------------------------------------------------------------------
n_states                    = 11;       % 10 base + 1 PDA flow state
X0                          = zeros(n_states, 1);

X0(idx.P_ra)    =  4;   % [mmHg] — RA mean pressure, neonate — Ref [3]
X0(idx.P_rv)    =  5;   % [mmHg] — RV diastolic, neonate      — Ref [3]
X0(idx.P_pa)    = 15;   % [mmHg] — PA mean, neonate (elevated) — Ref [3]
X0(idx.Q_pa_pul)=  0;   % [mL/s] — initial pulmonary flow
X0(idx.P_pv)    =  8;   % [mmHg] — pulmonary venous pressure   — Ref [3]
X0(idx.P_la)    =  5;   % [mmHg] — LA mean pressure, neonate   — Ref [3]
X0(idx.P_lv)    =  5;   % [mmHg] — LV diastolic pressure       — Ref [3]
X0(idx.P_ao)    = 50;   % [mmHg] — neonate aortic diastolic    — Ref [3]
X0(idx.Q_ao_sys)=  0;   % [mL/s] — initial systemic flow
X0(idx.P_sys)   = 50;   % [mmHg] — systemic venous pressure    — Ref [3]
X0(idx.Q_shunt_pda) = 0;% [mL/s] — initial PDA shunt flow

params.X0 = X0;

end
