function [E_lv, dE_lv, E_rv, dE_rv] = elastance_model(t, params)
% ELASTANCE_MODEL
% -----------------------------------------------------------------------
% Computes instantaneous left and right ventricular elastance and their
% time derivatives using the piecewise quadratic/cosine normalised model.
%
% The normalised activation function En(tn) follows:
%   Phase 1 (contraction,  0 ≤ tn ≤ Ts1): En = 1.55 * (tn/Ts1)²
%   Phase 2 (relaxation, Ts1 < tn ≤ Ts2): En = 1.55 * cos²(π/2*(tn-Ts1)/(Ts2-Ts1))
%   Phase 3 (diastole,   tn > Ts2):        En = 0
%
% Time-varying elastance:
%   E_lv(t) = (Emax_lv - Emin_lv) * En(tn) + Emin_lv
%   E_rv(t) = (Emax_rv - Emin_rv) * En(tn) + Emin_rv
%
% INPUTS:
%   t       - current simulation time                              [s]
%   params  - parameter struct (requires: T_cardiac, Ts1, Ts2,
%             Emax_lv, Emin_lv, Emax_rv, Emin_rv)
%
% OUTPUTS:
%   E_lv    - LV instantaneous elastance                          [mmHg/mL]
%   dE_lv   - dE_lv/dt                                           [mmHg/mL/s]
%   E_rv    - RV instantaneous elastance                          [mmHg/mL]
%   dE_rv   - dE_rv/dt                                           [mmHg/mL/s]
%
% ASSUMPTIONS:
%   - Both ventricles share the same normalised activation function En
%   - LV and RV contract synchronously (no inter-ventricular delay)
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%       Eqs. (5)–(7); piecewise approximation of Stergiopulos (1996)
%   [2] Stergiopulos N et al. (1996). Am J Physiol 270(6):H2050–H2059.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

% Normalised time within cardiac cycle [s]
tn = mod(t, params.T_cardiac);

%% Piecewise activation function En and its derivative dEn
if tn <= params.Ts1
    % Phase 1: contraction — quadratic ramp
    En  = 1.55 * (tn / params.Ts1)^2;
    dEn = 1.55 * 2 * (tn / params.Ts1) * (1 / params.Ts1);

elseif tn <= params.Ts2
    % Phase 2: relaxation — cosine decay
    alpha = (pi / 2) * (tn - params.Ts1) / (params.Ts2 - params.Ts1);
    En    = 1.55 * cos(alpha)^2;
    dEn   = 1.55 * 2 * cos(alpha) * (-sin(alpha)) * ...
            (pi / 2) / (params.Ts2 - params.Ts1);

else
    % Phase 3: diastole — fully relaxed
    En  = 0;
    dEn = 0;
end

%% Instantaneous LV elastance and derivative
E_lv  = (params.Emax_lv - params.Emin_lv) * En  + params.Emin_lv;  % [mmHg/mL]
dE_lv = (params.Emax_lv - params.Emin_lv) * dEn;                    % [mmHg/mL/s]

%% Instantaneous RV elastance and derivative
E_rv  = (params.Emax_rv - params.Emin_rv) * En  + params.Emin_rv;  % [mmHg/mL]
dE_rv = (params.Emax_rv - params.Emin_rv) * dEn;                    % [mmHg/mL/s]

end
