function [Q_tv, Q_pv_valve, Q_mv, Q_av] = valve_model(P_ra, P_rv, P_pa, P_la, P_lv, P_ao, params)
% VALVE_MODEL
% -----------------------------------------------------------------------
% Computes instantaneous flow through all four cardiac valves.
% Valves are modelled as ideal diodes with a fixed resistance, using
% a smooth max(0, ...) formulation to avoid integration discontinuities.
%
% Physics:
%   Q_valve = max(0, P_upstream - P_downstream) / R_valve
%
% INPUTS:
%   P_ra    - right atrial pressure                               [mmHg]
%   P_rv    - right ventricular pressure                          [mmHg]
%   P_pa    - pulmonary artery pressure                           [mmHg]
%   P_la    - left atrial pressure                                [mmHg]
%   P_lv    - left ventricular pressure                           [mmHg]
%   P_ao    - aortic pressure (proximal)                          [mmHg]
%   params  - parameter struct (R_tv, R_pv_valve, R_mv, R_av)
%
% OUTPUTS:
%   Q_tv        - tricuspid valve flow  (RA → RV, forward positive)  [mL/s]
%   Q_pv_valve  - pulmonary valve flow  (RV → PA, forward positive)  [mL/s]
%   Q_mv        - mitral valve flow     (LA → LV, forward positive)  [mL/s]
%   Q_av        - aortic valve flow     (LV → Ao, forward positive)  [mL/s]
%
% ASSUMPTIONS:
%   - Valve resistance is constant (no leaflet inertia model)
%   - max(0, dP)/R prevents backflow (ideal diode behaviour)
%   - No valve regurgitation or stenosis modelled at baseline
%
% SIGN CONVENTIONS:
%   Positive flow = physiologically forward direction for each valve.
%   Backflow is forced to zero by the diode (max) operator.
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%       Section 2.1 valve model.
%   [2] Guardrails §8.4 — smooth valve switching.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

% Tricuspid valve: RA → RV
Q_tv       = max(0, P_ra - P_rv) / params.R_tv;           % [mL/s]

% Pulmonary valve: RV → PA
Q_pv_valve = max(0, P_rv - P_pa) / params.R_pv_valve;     % [mL/s]

% Mitral valve: LA → LV
Q_mv       = max(0, P_la - P_lv) / params.R_mv;           % [mL/s]

% Aortic valve: LV → Aorta
Q_av       = max(0, P_lv - P_ao) / params.R_av;           % [mL/s]

end
