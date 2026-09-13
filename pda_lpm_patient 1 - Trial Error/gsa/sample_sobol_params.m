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
%     - Total evaluations: N × (D + 2)
%       (N for A, N for B, and N×D for AB_i)
%     - BA_i matrices are not generated.
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
    'R_pa',              0.05,     1.00   % [mmHg·s/mL] Pulmonary resistance — updated:
                                          %   clinical override (P_pa-P_pv)/Q_pul gives
                                          %   0.21–0.62 across patients; old [0.01,0.20]
                                          %   placed ALL samples below the operating point
    'C_ao',              0.0001,   0.0020 % [mL/mmHg]   Aortic compliance
    'C_sys',             0.05,     0.50   % [mL/mmHg]   Systemic compliance
    'Emax_lv',           3.0,      20.0   % [mmHg/mL]   LV elastance — updated:
                                          %   narrowed from [0.5,25] to neonatal clamp
                                          %   range [3,20]; wide range caused Emax_lv
                                          %   to monopolise variance (ST > 1 artefact)
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

%% 2. Generate reproducible Sobol base matrices
sampling_seed = 42;

rng_state_before_sampling = rng;
rng_cleanup = onCleanup(@() rng(rng_state_before_sampling));

rng(sampling_seed, 'twister');

try
    sob = sobolset(2 * D, 'Skip', 1, 'Leap', 31);
    sob = scramble(sob, 'MatousekAffineOwen');

    U = net(sob, N);

    U_A = U(:, 1:D);
    U_B = U(:, D+1:2*D);

    sampler_name = 'sobolset';
    scramble_method = 'MatousekAffineOwen';

catch ME
    warning('SAMPLE_SOBOL_PARAMS:SobolUnavailable', ...
        ['Sobol sampling unavailable (%s). ', ...
         'Falling back to pseudo-random sampling.'], ...
        ME.message);

    rng(sampling_seed, 'twister');

    U_A = rand(N, D);
    U_B = rand(N, D);

    sampler_name = 'rand';
    scramble_method = 'none';
end

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
sample_info.N                = N;
sample_info.D                = D;
sample_info.n_total          = n_total;
sample_info.idx_A            = (1:N)';
sample_info.idx_B            = (N+1:2*N)';
sample_info.idx_ABi          = idx_ABi;

% Reproducibility metadata
sample_info.sampling_seed    = sampling_seed;
sample_info.sampler_name     = sampler_name;
sample_info.scramble_method  = scramble_method;
sample_info.sobol_dimensions = 2 * D;

%% 7. Display parameter ranges
fprintf('\n  %-18s %12s %12s\n', 'Parameter', 'Lower', 'Upper');
fprintf('  %s\n', repmat('-', 1, 44));
for d = 1:D
    fprintf('  %-18s %12.4f %12.4f\n', param_names{d}, param_bounds(d,1), param_bounds(d,2));
end
fprintf('\n');

end
