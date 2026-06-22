function dX = system_rhs_pda_coa(t, X, params)
% SYSTEM_RHS_PDA_COA
% -----------------------------------------------------------------------
% ODE right-hand side for the 13-state cardiovascular LPM with PDA
% plus a virtual Coarctation of the Aorta (CoA) in the descending aorta.
%
% SCENARIO EXTENSION: Virtual CoA
%   State idx.P_ao_dist = 12: distal aortic pressure  [mmHg]
%   State idx.Q_coa     = 13: CoA segment flow        [mL/s]
%   This extension is active when params.scenario_coa == true.
%   Reference: docs/theory_notes.md, Section 1.3
%
% CoA PHYSICS:
%   The CoA segment is modelled as viscous (Poiseuille) + turbulent
%   (Bernoulli) resistance + inertance:
%     ΔP_coa = R_coa_viscous × Q_coa + K_turb_coa × Q_coa × |Q_coa|
%
%   CoA geometry parameters (R_coa_viscous, K_turb_coa, L_coa) are
%   computed in build_coa_params from:
%     - stenosis_pct  : anatomical area reduction           [%]
%     - coa_length_mm : CoA segment length (patient-specific or scenario)
%       Discrete/short CoA: < 5 mm  |  Long-segment CoA: ≥ 5 mm
%
% PRIMARY CLINICAL OUTPUT:
%   The predicted CoA severity is NOT taken from stenosis_pct directly.
%   It is classified in compute_clinical_indices from the simulated
%   mean-systolic ΔP_CoA against literature-based thresholds:
%     Mild     : ΔP_mean_sys <  20 mmHg
%     Moderate : ΔP_mean_sys 20–40 mmHg
%     Severe   : ΔP_mean_sys >  40 mmHg
%
% PDA ROLE:
%   PDA parameters are retained in params (R_shunt_pda, L_shunt_pda).
%   A patent PDA raises P_ao_dist, reducing observed ΔP_CoA and potentially
%   masking true CoA severity. Report PDA as a haemodynamic MODIFIER.
%
% State vector:
%   1  P_ra       [mmHg]  Right atrial pressure
%   2  P_rv       [mmHg]  Right ventricular pressure
%   3  P_pa       [mmHg]  Pulmonary artery pressure
%   4  Q_pa_pul   [mL/s]  Pulmonary arterial flow
%   5  P_pv       [mmHg]  Pulmonary venous pressure
%   6  P_la       [mmHg]  Left atrial pressure
%   7  P_lv       [mmHg]  Left ventricular pressure
%   8  P_ao       [mmHg]  Proximal aortic pressure
%   9  Q_ao_sys   [mL/s]  Systemic aortic flow (proximal segment)
%  10  P_sys      [mmHg]  Systemic venous pressure
%  11  Q_shunt_pda[mL/s]  PDA shunt flow (L→R positive)
%  12  P_ao_dist  [mmHg]  Distal aortic pressure (post-CoA)
%  13  Q_coa      [mL/s]  Flow through CoA segment
%
% INPUTS:
%   t       - current time                                         [s]
%   X       - state vector (13×1)
%   params  - parameter struct from build_coa_params
%             (must include coa_length_mm, coa_length_category,
%              severity_thresholds, R_coa_viscous, K_turb_coa, L_coa)
%
% OUTPUTS:
%   dX      - time derivatives (13×1)
%
% SIGN CONVENTIONS:
%   - Q_coa > 0       : forward flow, proximal → distal aorta
%   - Q_shunt_pda > 0 : left-to-right PDA shunt (Ao_dist → PA)
%   - DeltaP_coa      = P_ao (proximal) − P_ao_dist (distal)
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%   [2] Keshavarz-Motamed et al. (2011). J Biomech 44:2817–2825. Eq. 3–5.
%   [3] Baumgartner et al. (2010). Eur Heart J 31(19):2369–2417.
%       (ESC CoA severity gradient thresholds)
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  2.0  — variable CoA length; gradient-based severity
% -----------------------------------------------------------------------

dX  = zeros(13, 1);
idx = params.idx;

%% Unpack state vector
P_ra        = X(idx.P_ra);           % [mmHg]
P_rv        = X(idx.P_rv);           % [mmHg]
P_pa        = X(idx.P_pa);           % [mmHg]
Q_pa_pul    = X(idx.Q_pa_pul);       % [mL/s]
P_pv        = X(idx.P_pv);           % [mmHg]
P_la        = X(idx.P_la);           % [mmHg]
P_lv        = X(idx.P_lv);           % [mmHg]
P_ao        = X(idx.P_ao);           % [mmHg] — proximal aorta
Q_ao_sys    = X(idx.Q_ao_sys);       % [mL/s] — systemic (upper body) flow
P_sys       = X(idx.P_sys);          % [mmHg]
Q_shunt_pda = X(idx.Q_shunt_pda);    % [mL/s]
P_ao_dist   = X(idx.P_ao_dist);      % [mmHg] — distal aorta (post-CoA)
Q_coa       = X(idx.Q_coa);          % [mL/s] — CoA segment flow

%% 1. Ventricular elastance
[E_lv, dE_lv, E_rv, dE_rv] = elastance_model(t, params);

%% 2. Cardiac valve flows
[Q_tv, Q_pv_valve, Q_mv, Q_av] = valve_model(P_ra, P_rv, P_pa, P_la, P_lv, P_ao, params);

%% 3. Systemic flows
% Upper body (proximal, drains to systemic veins from proximal aorta)
Q_sys_upper = (P_ao - P_sys) / params.R_systemic;            % [mL/s]
% Lower body (distal, drains from distal aorta post-CoA)
Q_sys_lower = (P_ao_dist - P_sys) / params.R_systemic;       % [mL/s]
% Pulmonary venous return
Q_pv_return = (P_pv - P_la) / params.R_pv_veins;             % [mL/s]

%% 4. CoA pressure drop — viscous + turbulent (Ref [2] Eq. 3–5)
% DeltaP = R_viscous * Q + K_turb * Q * |Q|
DeltaP_coa_resistive = params.R_coa_viscous * Q_coa ...
                     + params.K_turb_coa * Q_coa * abs(Q_coa); % [mmHg]

%% 5. ODEs

% Right atrium
dX(idx.P_ra)       = (Q_sys_upper + Q_sys_lower - Q_tv) / params.C_ra;      % [mmHg/s]

% Right ventricle (active elastance)
dX(idx.P_rv)       = (P_rv / E_rv) * dE_rv + E_rv * (Q_tv - Q_pv_valve);   % [mmHg/s]

% Pulmonary artery: RV output + PDA inflow
dX(idx.P_pa)       = (Q_pv_valve + Q_shunt_pda - Q_pa_pul) / params.C_pa;  % [mmHg/s]

% Pulmonary flow momentum
dX(idx.Q_pa_pul)   = (P_pa - P_pv) / params.L_pa ...
                   - Q_pa_pul * params.R_pa / params.L_pa;                  % [mL/s²]

% Pulmonary veins
dX(idx.P_pv)       = (Q_pa_pul - Q_pv_return) / params.C_pv;               % [mmHg/s]

% Left atrium
dX(idx.P_la)       = (Q_pv_return - Q_mv) / params.C_la;                   % [mmHg/s]

% Left ventricle (active elastance)
dX(idx.P_lv)       = (P_lv / E_lv) * dE_lv + E_lv * (Q_mv - Q_av);       % [mmHg/s]

% Proximal aorta: LV output → feeds upper body + CoA segment (PDA no longer taps here)
dX(idx.P_ao)       = (Q_av - Q_ao_sys - Q_coa) / params.C_ao;               % [mmHg/s]

% Proximal systemic flow momentum
dX(idx.Q_ao_sys)   = (P_ao - P_sys) / params.L_ao ...
                   - Q_ao_sys * params.R_ao / params.L_ao;                  % [mL/s²]

% Systemic veins: receives both upper and lower body drainage
dX(idx.P_sys)      = (Q_ao_sys + Q_sys_lower - ...
                       (P_sys - P_ra) / params.R_systemic) / params.C_sys;  % [mmHg/s]

% PDA shunt flow momentum (Ao_dist → PA): PDA anatomically connects PA to descending aorta
dX(idx.Q_shunt_pda)= (P_ao_dist - P_pa) / params.L_shunt_pda ...
                   - Q_shunt_pda * params.R_shunt_pda / params.L_shunt_pda; % [mL/s²]

% Distal aorta: receives CoA flow, loses flow to PDA shunt and lower body
dX(idx.P_ao_dist)  = (Q_coa - Q_shunt_pda - Q_sys_lower) / params.C_sys;  % [mmHg/s]

% CoA flow momentum (proximal → distal through stenosis)
dX(idx.Q_coa)      = (P_ao - P_ao_dist - DeltaP_coa_resistive) / params.L_coa; % [mL/s²]

end
