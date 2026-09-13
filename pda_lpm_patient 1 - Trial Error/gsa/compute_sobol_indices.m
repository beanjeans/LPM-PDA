function [S1, ST, S1_ci, ST_ci] = compute_sobol_indices(Y, sample_info)
% COMPUTE_SOBOL_INDICES
% -----------------------------------------------------------------------
% Computes first-order (S1) and total-order (ST) Sobol sensitivity indices
% from model output Y evaluated on the Saltelli sampling scheme.
%
% METHOD:
%   Uses the Jansen (1999) / Saltelli (2010) estimators:
%
%   First-order (S1_i):
%     S1_i = [ Var(Y) - (1/2N) Σ (f(B)_j - f(AB_i)_j)² ] / Var(Y)
%          = 1 - V_i / Var(Y)
%     where V_i = (1/2N) Σ (f(B)_j - f(AB_i)_j)²
%
%   This is equivalent to the Saltelli (2010) Eq. (b), using the convention
%   that AB_i = A with column i taken from B (see sample_sobol_params.m):
%     V_i = (1/N) Σ f(B)_j × [f(AB_i)_j - f(A)_j]
%     S1_i = V_i / Var(Y)
%
%   Total-order (ST_i) — Jansen (1999) estimator:
%     VT_i = (1/2N) Σ (f(A)_j - f(AB_i)_j)²
%     ST_i = VT_i / Var(Y)
%
% CONFIDENCE INTERVALS:
%   Bootstrap resampling (n_boot=1000) for 95% confidence intervals.
%
% HANDLING NaN:
%   Rows where any output is NaN are excluded from variance computation.
%   If > 50% of samples are NaN, a warning is issued.
%
% INPUTS:
%   Y            - (n_total × n_outputs) matrix from evaluate_model_outputs
%   sample_info  - struct from sample_sobol_params containing:
%       .N, .D, .n_total, .idx_A, .idx_B, .idx_ABi
%
% OUTPUTS:
%   S1     - (D × n_outputs) first-order Sobol indices
%   ST     - (D × n_outputs) total-order Sobol indices
%   S1_ci  - (D × n_outputs × 2) confidence intervals [lower, upper]
%   ST_ci  - (D × n_outputs × 2) confidence intervals [lower, upper]
%
% REFERENCES:
%   [1] Saltelli A et al. (2010). Variance based sensitivity analysis of
%       model output. Comp Phys Comm 181:259–270. Eqs. (b) and (f).
%   [2] Jansen MJW (1999). Analysis of variance designs for model output.
%       Comp Phys Comm 117:35–43.
%   [3] Sobol IM (2001). Global sensitivity indices for nonlinear
%       mathematical models. Math Comp Simul 55:271–280.
%
% AUTHOR:   GSA Extension — Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

N = sample_info.N;
D = sample_info.D;
n_outputs = size(Y, 2);

%% 1. Extract output vectors for each block
% -----------------------------------------------------------------------
Y_A = Y(sample_info.idx_A, :);     % N × n_outputs
Y_B = Y(sample_info.idx_B, :);     % N × n_outputs

% Y_ABi{i} = N × n_outputs for the AB_i matrix
Y_ABi = cell(1, D);
for i = 1:D
    Y_ABi{i} = Y(sample_info.idx_ABi(:, i), :);
end

%% 2. Compute indices for each output
S1 = zeros(D, n_outputs);
ST = zeros(D, n_outputs);

n_boot = 1000;
S1_ci = zeros(D, n_outputs, 2);
ST_ci = zeros(D, n_outputs, 2);

n_valid_per_output = zeros(1, n_outputs);  % Track valid sample count per output

bootstrap_seed = 12345;

% Dedicated bootstrap stream. This does not modify MATLAB's global RNG.
bootstrap_stream = RandStream('mt19937ar', ...
    'Seed', bootstrap_seed);

for q = 1:n_outputs

    y_A = Y_A(:, q);
    y_B = Y_B(:, q);

    % Collect AB_i columns for this output
    y_ABi_q = zeros(N, D);
    for i = 1:D
        y_ABi_q(:, i) = Y_ABi{i}(:, q);
    end

    % --- Handle NaN: find rows valid across A, B, and all AB_i -----------
    valid_mask = ~isnan(y_A) & ~isnan(y_B);
    for i = 1:D
        valid_mask = valid_mask & ~isnan(y_ABi_q(:, i));
    end
    n_valid = sum(valid_mask);
    n_valid_per_output(q) = n_valid;   % Record for summary

    if n_valid < 0.5 * N
        warning('COMPUTE_SOBOL_INDICES: Output %d has %d/%d valid samples (<50%%). Results may be unreliable.', ...
            q, n_valid, N);
    end

    if n_valid < 10
        warning('COMPUTE_SOBOL_INDICES: Output %d has too few valid samples (%d). Skipping.', q, n_valid);
        S1(:, q) = NaN;
        ST(:, q) = NaN;
        S1_ci(:, q, :) = NaN;
        ST_ci(:, q, :) = NaN;
        continue;
    end

    % Apply mask
    yA  = y_A(valid_mask);
    yB  = y_B(valid_mask);
    yAB = y_ABi_q(valid_mask, :);
    Nv  = n_valid;

    % --- Total variance (pooled from A and B) ---
    y_all_AB = [yA; yB];
    f0  = mean(y_all_AB);          % Grand mean
    VarY = var(y_all_AB, 1);       % Population variance (divide by n)

    if VarY < eps
        % Constant output → all indices are zero
        S1(:, q) = 0;
        ST(:, q) = 0;
        continue;
    end

    % --- Compute S1 and ST for each parameter ---
    for i = 1:D
        yABi = yAB(:, i);

        % --- First-order index (Saltelli 2010 estimator, Eq. b) ---
        % AB_i = A with column i taken from B (see sample_sobol_params.m).
        % For this convention the correct pairing is:
        %   V_i = (1/N) * sum( f(B) * (f(AB_i) - f(A)) )
        %   S1_i = V_i / Var(Y)
        V_first = (1/Nv) * sum(yB .* (yABi - yA));
        S1(i, q) = V_first / VarY;

        % --- Total-order index (Jansen 1999 estimator) ---
        % VT_i = (1/2N) * sum( (f(A) - f(AB_i))^2 )
        V_total = (1/(2*Nv)) * sum((yA - yABi).^2);
        ST(i, q) = V_total / VarY;
    end

    % --- Bootstrap confidence intervals ---
    S1_boot = zeros(D, n_boot);
    ST_boot = zeros(D, n_boot);

    for b = 1:n_boot
        % Resample indices with replacement
        boot_idx = randi(bootstrap_stream, Nv, Nv, 1);
        yA_b  = yA(boot_idx);
        yB_b  = yB(boot_idx);

        y_all_b = [yA_b; yB_b];
        VarY_b  = var(y_all_b, 1);

        if VarY_b < eps
            S1_boot(:, b) = 0;
            ST_boot(:, b) = 0;
            continue;
        end

        for i = 1:D
            yABi_b = yAB(boot_idx, i);

            V_first_b = (1/Nv) * sum(yB_b .* (yABi_b - yA_b));
            S1_boot(i, b) = V_first_b / VarY_b;

            V_total_b = (1/(2*Nv)) * sum((yA_b - yABi_b).^2);
            ST_boot(i, b) = V_total_b / VarY_b;
        end
    end

    % 95% CI from bootstrap percentiles
    for i = 1:D
        S1_ci(i, q, 1) = prctile(S1_boot(i, :), 2.5);
        S1_ci(i, q, 2) = prctile(S1_boot(i, :), 97.5);
        ST_ci(i, q, 1) = prctile(ST_boot(i, :), 2.5);
        ST_ci(i, q, 2) = prctile(ST_boot(i, :), 97.5);
    end
end

%% 3. Display summary
fprintf('\n=== SOBOL INDICES COMPUTED ===\n');
fprintf('  Valid samples (per output): min=%d  max=%d  (of N=%d total)\n', ...
    min(n_valid_per_output), max(n_valid_per_output), N);
fprintf('  Bootstrap samples:   %d (95%% CI)\n', n_boot);
fprintf('\n');

end
