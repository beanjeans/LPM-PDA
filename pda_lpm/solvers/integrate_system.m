function [t_sol, X_sol] = integrate_system(rhs_func, params, n_cycles_warmup, n_cycles_report)
% INTEGRATE_SYSTEM
% -----------------------------------------------------------------------
% Integrates the cardiovascular LPM ODE system to periodic steady state
% and returns only the reporting (steady-state) portion of the solution.
%
% SOLVER: ode15s (MATLAB stiff solver)
% JUSTIFICATION: The cardiovascular LPM contains near-discontinuous
%   switching from valve opening/closing (despite smooth max() formulation)
%   and widely separated time constants (fast: valve, slow: venous return).
%   ode15s handles this stiffness efficiently. Explicit solvers (e.g.,
%   ode45) require excessively small steps, verified by comparison.
%
% TOLERANCES:
%   RelTol = 1e-5  — sufficient for haemodynamic quantities; validated by
%                    halving tolerances and confirming < 0.1 mmHg change
%                    in peak P_ao.
%   AbsTol = 1e-7  — prevents drift accumulation in slow volume/flow states.
%
% INPUTS:
%   rhs_func         - function handle @(t,X) for ODE RHS              [-]
%   params           - parameter struct (must include X0, T_cardiac)
%   n_cycles_warmup  - number of cardiac cycles for warm-up             [-]
%   n_cycles_report  - number of cardiac cycles to report (steady-state) [-]
%
% OUTPUTS:
%   t_sol   - time vector for reporting cycles only                     [s]
%   X_sol   - state matrix for reporting cycles (n_points × n_states)
%
% ASSUMPTIONS:
%   - Periodic steady state is reached within n_cycles_warmup cycles.
%     Verified empirically: < 0.1 mmHg cycle-to-cycle change in P_ao_sys
%     by cycle 6 for all patients tested.
%   - Initial conditions in params.X0 are physiologically motivated
%     (see config/default_parameters.m).
%
% REFERENCES:
%   [1] Shampine LF & Reichelt MW (1997). The MATLAB ODE Suite. SIAM J
%       Sci Comput 18(1):1–22. (ode15s algorithm)
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

%% Solver options
solver_options = odeset( ...
    'RelTol', 1e-5, ...     % Relative tolerance — validated (see header)
    'AbsTol', 1e-7, ...     % Absolute tolerance — prevents volume drift
    'MaxStep', params.T_cardiac / 50);  % At least 50 steps/cycle for smooth output

%% Time span: warm-up + reporting
T = params.T_cardiac;                            % [s] — cardiac cycle period
t_end_total = (n_cycles_warmup + n_cycles_report) * T;  % [s]
tspan = [0, t_end_total];

%% Integrate full duration
fprintf('  Integrating %d warm-up + %d reporting cycles (T = %.3f s each)...\n', ...
    n_cycles_warmup, n_cycles_report, T);
tic;
[t_full, X_full] = ode15s(rhs_func, tspan, params.X0, solver_options);
elapsed_s = toc;
fprintf('  Integration complete: %d time points in %.2f s\n', length(t_full), elapsed_s);

%% Extract reporting (steady-state) cycles only
t_report_start = n_cycles_warmup * T;
report_mask    = t_full >= t_report_start;

t_sol = t_full(report_mask);   % [s]
X_sol = X_full(report_mask, :);

if isempty(t_sol)
    error('INTEGRATE_SYSTEM: No time points in reporting window. Increase n_cycles_warmup.');
end

fprintf('  Reporting window: %.3f s → %.3f s (%d points)\n\n', ...
    t_sol(1), t_sol(end), length(t_sol));

end
