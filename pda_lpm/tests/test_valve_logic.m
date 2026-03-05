%% TEST_VALVE_LOGIC
% -----------------------------------------------------------------------
% Unit tests for valve_model.m — verifies that all valve flow directions
% and zero-backflow conditions are correct.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

clear; clc;
addpath('../config', '../models', '../utils');

fprintf('=================================================================\n');
fprintf('  TEST_VALVE_LOGIC\n');
fprintf('=================================================================\n\n');

params = default_parameters();
n_pass = 0; n_fail = 0;

function result = assert_approx(name, val, expected, tol)
    err = abs(val - expected);
    if err <= tol
        fprintf('  ✓ %s = %.4f (expected %.4f)\n', name, val, expected);
        result = 1;
    else
        fprintf('  ✗ %s = %.4f (expected %.4f, error = %.4f)\n', name, val, expected, err);
        result = 0;
    end
end

%% Test 1: Forward flow through all valves
[Q_tv, Q_pv, Q_mv, Q_av] = valve_model(10, 5, 3, 8, 4, 70, params);
n_pass = n_pass + assert_approx('Q_tv (forward)', Q_tv, (10-5)/params.R_tv, 1e-10);
n_pass = n_pass + assert_approx('Q_pv (forward)', Q_pv, (5-3)/params.R_pv_valve, 1e-10);
n_pass = n_pass + assert_approx('Q_mv (forward)', Q_mv, (8-4)/params.R_mv, 1e-10);
n_pass = n_pass + assert_approx('Q_av (forward)', Q_av, (4-70)/params.R_av, 1e-10);

%% Test 2: Zero backflow — pressure reversal
[Q_tv2, Q_pv2, Q_mv2, Q_av2] = valve_model(3, 10, 20, 3, 10, 80, params);
n_pass = n_pass + assert_approx('Q_tv (closed, no backflow)', Q_tv2, 0, 1e-10);
n_pass = n_pass + assert_approx('Q_pv (closed, no backflow)', Q_pv2, 0, 1e-10);
n_pass = n_pass + assert_approx('Q_mv (closed, no backflow)', Q_mv2, 0, 1e-10);
n_pass = n_pass + assert_approx('Q_av (closed, no backflow)', Q_av2, 0, 1e-10);

fprintf('\n--- Summary ---\n');
fprintf('  Passed: %d / %d\n', n_pass, n_pass + n_fail);
if n_fail == 0
    fprintf('  STATUS: ALL VALVE LOGIC TESTS PASSED\n');
else
    fprintf('  STATUS: %d FAILED\n', n_fail);
end
fprintf('=================================================================\n\n');
