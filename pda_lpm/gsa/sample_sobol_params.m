function [X_all, param_names, param_bounds, n_total, sample_info] = sample_sobol_params(N)
% SAMPLE_SOBOL_PARAMS
% -----------------------------------------------------------------------
% Generates Saltelli sampling matrices for Sobol Global Sensitivity
% Analysis of the PDA-CoA LPM model.
%
% METHOD:
%   Uses the Saltelli (2002) sampling scheme:
%     - Generate two independent quasi-random base matrices A, B of
%       size N × D in [0,1]^D using Sobol quasi-random sequences.
%     - For each parameter i (1..D), construct AB_i by taking column i
%       from B and all other columns from A.
%     - Total evaluations: N × (2D + 2)
%       (N for A, N for B, N×D for AB_i, N×D for BA_i)
%       Simplified: we use A, B, and D copies of AB_i → N*(D+2)
%
% PARAMETER SET (reduced, 8 parameters):
%   1. R_shunt_pda    — PDA shunt resistance        [mmHg·s/mL]
%   2. R_systemic     — Total systemic resistance    [mmHg·s/mL]
%   3. R_pa           — Pulmonary arterial resistance[mmHg·s/mL]
%   4. C_ao           — Aortic compliance            [mL/mmHg]
%   5. C_sys          — Systemic compliance           [mL/mmHg]
%   6. Emax_lv        — LV peak elastance            [mmHg/mL]
%   7. stenosis_pct   — CoA stenosis severity        [%]
%   8. coa_length_mm  — CoA segment length           [mm]
%
% INPUTS:
%   N  - base sample size (e.g. 256 or 512)
%        Total model evaluations = N × (D + 2) where D = 8
%        → N=256: 2560 evaluations
%        → N=512: 5120 evaluations
%
% OUTPUTS:
%   X_all        - (n_total × D) matrix of parameter samples,
%                  each row is one parameter combination in physical units
%   param_names  - (1 × D) cell array of parameter names
%   param_bounds - (D × 2) matrix of [lower, upper] bounds
%   n_total      - total number of model evaluations
%   sample_info  - struct with metadata:
%       .N         — base sample size
%       .D         — number of parameters
%       .n_total   — total evaluations
%       .idx_A     — row indices for matrix A evaluations
%       .idx_B     — row indices for matrix B evaluations
%       .idx_ABi   — (N × D) matrix: idx_ABi(:,i) = rows for AB_i
%
% REFERENCES:
%   [1] Saltelli A (2002). Making best use of model evaluations to
%       compute sensitivity indices. Comp Phys Comm 145:280–297.
%   [2] Sobol IM (2001). Global sensitivity indices for nonlinear
%       mathematical models and their Monte Carlo estimates.
%       Math Comp Simul 55:271–280.
%   [3] Saltelli A et al. (2010). Variance based sensitivity analysis
%       of model output. Comp Phys Comm 181:259–270.
%
% AUTHOR:   GSA Extension — Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

%% 1. Define parameter names and bounds
% -----------------------------------------------------------------------
% Bounds chosen from physiological literature and neonatal ranges.
% These define the uncertainty/variability space for GSA.
%
%                         Name              Lower    Upper     Units
param_defs = {
    'R_shunt_pda',       0.5,      20.0   % [mmHg·s/mL] PDA resistance
    'R_systemic',        1.0,      15.0   % [mmHg·s/mL] Total SVR
    'R_pa',              0.01,     0.20   % [mmHg·s/mL] Pulmonary resistance
    'C_ao',              0.005,    0.05   % [mL/mmHg]   Aortic compliance
    'C_sys',             0.05,     0.50   % [mL/mmHg]   Systemic compliance
    'Emax_lv',           0.5,      5.0    % [mmHg/mL]   LV elastance
    'stenosis_pct',      10,       95     % [%]         CoA stenosis
    'coa_length_mm',     1,        15     % [mm]        CoA length
};

param_names  = param_defs(:, 1)';           % 1×D cell
param_bounds = cell2mat(param_defs(:, 2:3)); % D×2
D = size(param_bounds, 1);

fprintf('=== SOBOL SAMPLING ===\n');
fprintf('  Parameters (D):     %d\n', D);
fprintf('  Base sample (N):    %d\n', N);
fprintf('  Total evaluations:  N × (D+2) = %d\n', N * (D + 2));

%% 2. Generate quasi-random base samples in [0,1]^D
% -----------------------------------------------------------------------
% Use MATLAB's Sobol quasi-random sequence (sobolset) for low-discrepancy
% sampling. Skip the first point (origin) and leap for better uniformity.
% -----------------------------------------------------------------------
try
    sob = sobolset(D, 'Skip', 1, 'Leap', 31);
    sob = scramble(sob, 'MatousekAffineOwen');  % Randomise for unbiased estimation
    U = net(sob, 2 * N);  % Draw 2N points → split into A and B
    fprintf('  Sampler: sobolset (quasi-random, scrambled)\n');
catch
    % Fallback: if Statistics Toolbox not available, use pseudo-random
    warning('SAMPLE_SOBOL_PARAMS: sobolset unavailable — falling back to rand().');
    rng(42, 'twister');  % Fixed seed for reproducibility
    U = rand(2 * N, D);
    fprintf('  Sampler: rand (pseudo-random, seed=42)\n');
end

%% 3. Split into base matrices A and B
U_A = U(1:N, :);       % N × D in [0,1]
U_B = U(N+1:2*N, :);   % N × D in [0,1]

%% 4. Construct Saltelli sampling scheme
% -----------------------------------------------------------------------
% Stack: [A; B; AB_1; AB_2; ... AB_D]
% Where AB_i = A with column i replaced by B(:,i)
% Total rows: N + N + N*D = N*(D+2)
% -----------------------------------------------------------------------
n_total = N * (D + 2);

% Allocate in unit hypercube
U_all = zeros(n_total, D);

% Block 1: Matrix A
U_all(1:N, :) = U_A;

% Block 2: Matrix B
U_all(N+1:2*N, :) = U_B;

% Block 3..D+2: AB_i matrices
idx_ABi = zeros(N, D);  % Store row indices for each AB_i
for i = 1:D
    row_start = 2*N + (i-1)*N + 1;
    row_end   = 2*N + i*N;

    AB_i = U_A;               % Start from A
    AB_i(:, i) = U_B(:, i);   % Replace column i with B's column i

    U_all(row_start:row_end, :) = AB_i;
    idx_ABi(:, i) = (row_start:row_end)';
end

%% 5. Scale from [0,1] to physical parameter bounds
% -----------------------------------------------------------------------
% X_phys = lower + U * (upper - lower)
% -----------------------------------------------------------------------
X_all = zeros(n_total, D);
for d = 1:D
    lb = param_bounds(d, 1);
    ub = param_bounds(d, 2);
    X_all(:, d) = lb + U_all(:, d) * (ub - lb);
end

%% 6. Pack metadata
sample_info.N       = N;
sample_info.D       = D;
sample_info.n_total = n_total;
sample_info.idx_A   = (1:N)';
sample_info.idx_B   = (N+1:2*N)';
sample_info.idx_ABi = idx_ABi;

%% 7. Display parameter ranges
fprintf('\n  %-18s %12s %12s\n', 'Parameter', 'Lower', 'Upper');
fprintf('  %s\n', repmat('-', 1, 44));
for d = 1:D
    fprintf('  %-18s %12.4f %12.4f\n', param_names{d}, param_bounds(d,1), param_bounds(d,2));
end
fprintf('\n');

end
