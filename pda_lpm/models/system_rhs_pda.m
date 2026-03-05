function dX = system_rhs_pda(t, X, params)
% SYSTEM_RHS_PDA
% -----------------------------------------------------------------------
% ODE right-hand side for the 11-state cardiovascular LPM with PDA.
% Implements the Windkessel circuit of Ortiz-Rangel et al. (2022) plus
% a PDA shunt connecting the aorta to the pulmonary artery.
%
% State vector (access via params.idx — never hardcode indices):
%   X(idx.P_ra)        [mmHg]  Right atrial pressure
%   X(idx.P_rv)        [mmHg]  Right ventricular pressure
%   X(idx.P_pa)        [mmHg]  Pulmonary artery pressure
%   X(idx.Q_pa_pul)    [mL/s]  Pulmonary arterial flow
%   X(idx.P_pv)        [mmHg]  Pulmonary venous pressure
%   X(idx.P_la)        [mmHg]  Left atrial pressure
%   X(idx.P_lv)        [mmHg]  Left ventricular pressure
%   X(idx.P_ao)        [mmHg]  Aortic (proximal) pressure
%   X(idx.Q_ao_sys)    [mL/s]  Systemic aortic flow
%   X(idx.P_sys)       [mmHg]  Systemic venous pressure
%   X(idx.Q_shunt_pda) [mL/s]  PDA shunt flow (L→R positive)
%
% INPUTS:
%   t       - current time                                         [s]
%   X       - state vector (11×1)                                  [mixed]
%   params  - parameter struct from build_patient_params
%
% OUTPUTS:
%   dX      - time derivative of state vector (11×1)               [mixed/s]
%
% SIGN CONVENTIONS:
%   - Q_shunt_pda > 0 : flow from aorta to PA (left-to-right)
%   - All valve flows: positive = forward physiological direction
%   - Pressure gradients: always upstream minus downstream
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%       Eqs. (8)–(20).
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

dX = zeros(11, 1);
idx = params.idx;

%% Unpack state vector via index struct (Guardrail §7.1)
P_ra        = X(idx.P_ra);        % [mmHg]
P_rv        = X(idx.P_rv);        % [mmHg]
P_pa        = X(idx.P_pa);        % [mmHg]
Q_pa_pul    = X(idx.Q_pa_pul);    % [mL/s]
P_pv        = X(idx.P_pv);        % [mmHg]
P_la        = X(idx.P_la);        % [mmHg]
P_lv        = X(idx.P_lv);        % [mmHg]
P_ao        = X(idx.P_ao);        % [mmHg]
Q_ao_sys    = X(idx.Q_ao_sys);    % [mL/s]
P_sys       = X(idx.P_sys);       % [mmHg]
Q_shunt_pda = X(idx.Q_shunt_pda); % [mL/s]

%% 1. Time-varying ventricular elastance
[E_lv, dE_lv, E_rv, dE_rv] = elastance_model(t, params);

%% 2. Cardiac valve flows (smooth diode model)
[Q_tv, Q_pv_valve, Q_mv, Q_av] = valve_model(P_ra, P_rv, P_pa, P_la, P_lv, P_ao, params);

%% 3. Vascular bed flows
% Systemic venous return: from systemic veins to right atrium
Q_sys_return = (P_sys - P_ra) / params.R_systemic;       % [mL/s]

% Pulmonary venous return: from pulmonary veins to left atrium
Q_pv_return  = (P_pv - P_la) / params.R_pv_veins;        % [mL/s]

%% 4. PDA shunt flow contribution
% Q_shunt_pda is integrated as a state (momentum eq. below).
% Positive = aorta → PA (left-to-right, physiological in PDA with P_ao > P_pa)

%% 5. ODEs — Ortiz-Rangel (2022) Eqs. (8)–(17) + PDA Eqs. (18)–(20)

% Right atrium: receives systemic venous return, loses tricuspid flow
dX(idx.P_ra)     = (Q_sys_return - Q_tv) / params.C_ra;         % [mmHg/s]

% Right ventricle: active elastance formulation (Eq. 9)
dX(idx.P_rv)     = (P_rv / E_rv) * dE_rv ...
                 + E_rv * (Q_tv - Q_pv_valve);                   % [mmHg/s]

% Pulmonary artery: receives RV output + PDA inflow (Eq. 18)
dX(idx.P_pa)     = (Q_pv_valve + Q_shunt_pda - Q_pa_pul) / params.C_pa; % [mmHg/s]

% Pulmonary flow (inertial momentum, Eq. 11)
dX(idx.Q_pa_pul) = (P_pa - P_pv) / params.L_pa ...
                 - Q_pa_pul * params.R_pa / params.L_pa;         % [mL/s²]

% Pulmonary veins
dX(idx.P_pv)     = (Q_pa_pul - Q_pv_return) / params.C_pv;      % [mmHg/s]

% Left atrium
dX(idx.P_la)     = (Q_pv_return - Q_mv) / params.C_la;          % [mmHg/s]

% Left ventricle: active elastance formulation (Eq. 14)
dX(idx.P_lv)     = (P_lv / E_lv) * dE_lv ...
                 + E_lv * (Q_mv - Q_av);                         % [mmHg/s]

% Aorta: receives LV output, loses systemic flow + PDA shunt (Eq. 19)
dX(idx.P_ao)     = (Q_av - Q_ao_sys - Q_shunt_pda) / params.C_ao; % [mmHg/s]

% Systemic flow (inertial momentum, Eq. 16)
dX(idx.Q_ao_sys) = (P_ao - P_sys) / params.L_ao ...
                 - Q_ao_sys * params.R_ao / params.L_ao;          % [mL/s²]

% Systemic veins
dX(idx.P_sys)    = (Q_ao_sys - Q_sys_return) / params.C_sys;    % [mmHg/s]

% PDA shunt flow momentum (Eq. 20): positive = aorta → PA
dX(idx.Q_shunt_pda) = (P_ao - P_pa) / params.L_shunt_pda ...
                    - Q_shunt_pda * params.R_shunt_pda / params.L_shunt_pda; % [mL/s²]

end
