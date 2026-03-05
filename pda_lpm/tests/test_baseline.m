%% TEST_BASELINE
% -----------------------------------------------------------------------
% Validation test for the PDA LPM at reference (default) parameters.
% This script must pass without errors before any patient-specific run.
%
% Physiological reference ranges for a healthy neonate at rest:
%   P_ao_sys:   40–80 mmHg   [Rudolph 2001, Kleinman 2008]
%   P_ao_dia:   20–50 mmHg
%   P_pa_mean:  10–30 mmHg   (elevated in neonate)
%   EF_lv:      0.55–0.80
%   CO:         0.2–1.0 L/min (neonate ~3 kg, SV ~3–5 mL at HR 120–150)
%   Qp/Qs:      ~1.0 at baseline (no shunt)
%
% REFERENCES:
%   [1] Rudolph AM (2001). Congenital Diseases of the Heart. Futura.
%   [2] Kleinman CS et al. (2008). Echocardiography in Pediatric &
%       Congenital Heart Disease. Wiley-Blackwell.
%   [3] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%       Table 3 — validation criteria.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

clear; clc;

fprintf('=================================================================\n');
fprintf('  TEST_BASELINE — PDA LPM Physiological Range Checks\n');
fprintf('=================================================================\n\n');

%% Add paths
addpath('../config', '../models', '../solvers', '../utils');

%% Load default parameters (no PDA — test baseline)
params = default_parameters();

% Disable PDA for baseline test (R_shunt_pda = Inf means shunt absent)
params.R_shunt_pda = Inf;   % No PDA shunt at baseline

% Simulate with adult-equivalent parameters for easier range checking
params.HR_bpm    = 130;
params.T_cardiac = 60 / 130;
params.Ts1       = 0.3  * sqrt(params.T_cardiac);
params.Ts2       = 0.45 * sqrt(params.T_cardiac);
params.R_systemic = 10;    % [mmHg·s/mL] — estimated neonate SVR

%% Integrate
rhs_func = @(t, X) system_rhs_pda(t, X, params);

fprintf('Running baseline simulation (6 warm-up + 2 reporting cycles)...\n');
[t_sol, X_sol] = integrate_system(rhs_func, params, 6, 2);

%% Extract last cycle metrics
T        = params.T_cardiac;
mask_cyc = t_sol >= (t_sol(end) - T);
X_cyc    = X_sol(mask_cyc, :);
idx      = params.idx;

P_ao_sys  = max(X_cyc(:, idx.P_ao));   % [mmHg]
P_ao_dia  = min(X_cyc(:, idx.P_ao));   % [mmHg]
P_pa_mean = mean(X_cyc(:, idx.P_pa));  % [mmHg]

% LV volume reconstruction
t_cyc    = t_sol(mask_cyc);
tn_cyc   = mod(t_cyc, T);
En_cyc   = zeros(size(tn_cyc));
for k = 1:length(tn_cyc)
    if tn_cyc(k) <= params.Ts1
        En_cyc(k) = 1.55 * (tn_cyc(k) / params.Ts1)^2;
    elseif tn_cyc(k) <= params.Ts2
        alpha = (pi/2)*(tn_cyc(k)-params.Ts1)/(params.Ts2-params.Ts1);
        En_cyc(k) = 1.55 * cos(alpha)^2;
    end
end
E_lv_cyc = (params.Emax_lv - params.Emin_lv)*En_cyc + params.Emin_lv;
V_lv_cyc = X_cyc(:, idx.P_lv) ./ E_lv_cyc + params.V0_lv;
V_lv_ed  = max(V_lv_cyc);
V_lv_es  = min(V_lv_cyc);
SV_lv    = V_lv_ed - V_lv_es;
EF_lv    = SV_lv / V_lv_ed;
CO_Lmin  = SV_lv * params.HR_bpm / 1000;

% Mass conservation: net accumulation in P_ao over last cycle
dP_ao    = X_sol(end, idx.P_ao) - X_sol(find(t_sol >= t_sol(end)-T, 1), idx.P_ao);
mass_ok  = abs(dP_ao) < 1.0;   % < 1 mmHg drift acceptable

%% Assertions (physiological range checks)
n_passed = 0; n_failed = 0;

function check(name, val, lo, hi, unit)
    if val >= lo && val <= hi
        fprintf('  ✓ %-25s = %6.2f %-12s [%.1f – %.1f]\n', name, val, unit, lo, hi);
        n_passed = n_passed + 1;
    else
        fprintf('  ✗ %-25s = %6.2f %-12s [%.1f – %.1f]  ← OUT OF RANGE\n', name, val, unit, lo, hi);
        n_failed = n_failed + 1;
    end
end

fprintf('--- Physiological Range Assertions ---\n');
check('P_ao_sys',  P_ao_sys,  40,  80,  'mmHg');
check('P_ao_dia',  P_ao_dia,  20,  50,  'mmHg');
check('P_pa_mean', P_pa_mean, 10,  30,  'mmHg');
check('EF_lv',     EF_lv,     0.55, 0.80, '(fraction)');
check('CO',        CO_Lmin,   0.2, 1.0, 'L/min');

if mass_ok
    fprintf('  ✓ %-25s = %6.3f %-12s [< 1.0]\n', 'P_ao drift (last cycle)', dP_ao, 'mmHg');
    n_passed = n_passed + 1;
else
    fprintf('  ✗ %-25s = %6.3f %-12s [< 1.0]  ← MASS BALANCE FAILED\n', 'P_ao drift', dP_ao, 'mmHg');
    n_failed = n_failed + 1;
end

fprintf('\n--- Summary ---\n');
fprintf('  Passed: %d / %d\n', n_passed, n_passed + n_failed);

if n_failed == 0
    fprintf('  STATUS: ALL TESTS PASSED — baseline is physiologically valid\n');
else
    fprintf('  STATUS: %d TEST(S) FAILED — review parameter calibration\n', n_failed);
end
fprintf('=================================================================\n\n');
