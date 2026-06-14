%% RUN_SOBOL_GSA_PDA_COA
% =========================================================================
% GLOBAL SENSITIVITY ANALYSIS — SOBOL METHOD
% PDA-CoA Lumped Parameter Cardiovascular Model
%
% PURPOSE:
%   Perform Sobol/Saltelli Global Sensitivity Analysis (GSA) to identify
%   which model parameters most strongly influence the clinical CoA outputs:
%     1. ΔP_CoA_peak          — peak CoA pressure gradient         [mmHg]
%     2. ΔP_CoA_mean_systolic — mean systolic gradient             [mmHg]
%     3. Q_coa / Q_total      — CoA flow fra
% 11
% ction                  [0–1]
%     4. predicted_CoA_severity — clinical classification (1/2/3)
%
% WORKFLOW:
%   1. Generate Saltelli sampling matrices (quasi-random Sobol sequences)
%   2. Evaluate PDA-CoA model for each parameter combination
%   3. Compute first-order (S1) and total-order (ST) Sobol indices
%   4. Rank parameters by importance
%   5. Identify influential parameters for L-BFGS-B optimisation
%   6. Fix low-sensitivity parameters to default/literature values
%
% PARAMETERS ANALYSED (D=8):
%   R_shunt_pda, R_systemic, R_pa, C_ao, C_sys, Emax_lv,
%   stenosis_pct, coa_length_mm
%
% SAMPLING:
%   N = 256 (default) or 512
%   Total evaluations: N × (D + 2) = 2560 or 5120
%
% OUTPUT FILES (saved to pda_lpm/results/gsa/):
%   sobol_indices.csv          — S1, ST, and CIs for all outputs
%   sobol_ranking.csv          — parameters ranked by total-order index
%   gsa_raw_data.mat           — full workspace for reproducibility
%   sobol_barplot_S1.png       — bar chart of first-order indices
%   sobol_barplot_ST.png       — bar chart of total-order indices
%   sobol_barplot_combined.png — grouped S1+ST comparison
%
% NOTES:
%   - Does NOT modify any existing model files.
%   - Calls existing functions: build_coa_params, system_rhs_pda_coa,
%     integrate_system, compute_clinical_indices via evaluate_model_outputs.
%   - Patient data loaded non-interactively (batch mode).
%
% REFERENCES:
%   [1] Saltelli A (2002). Comp Phys Comm 145:280–297.
%   [2] Sobol IM (2001). Math Comp Simul 55:271–280.
%   [3] Saltelli A et al. (2010). Comp Phys Comm 181:259–270.
%
% AUTHOR:   GSA Extension — Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% =========================================================================

clear; clc; close all;

%% Add paths (existing code + new GSA directory)
addpath('config', 'models', 'solvers', 'utils', 'tests', 'gsa');

fprintf('=================================================================\n');
fprintf('   SOBOL GLOBAL SENSITIVITY ANALYSIS\n');
fprintf('   PDA-CoA Lumped Parameter Cardiovascular Model\n');
fprintf('=================================================================\n\n');

%% =========================================================================
%  CONFIGURATION
% =========================================================================
% Base sample size (N). Total evaluations = N * (D + 2) where D = 8.
%   N = 256 → 2560 evaluations (~30-45 min)
%   N = 512 → 5120 evaluations (~60-90 min)
N = 1024;

% Patient index (1-based, corresponds to patient_data.csv row)
% Set to desired patient or loop over all patients
patient_idx = 1;

% CSV path
csv_path = fullfile('config', 'patient_data.csv');

% Output directory
results_dir = fullfile('results', 'gsa');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('  Created output directory: %s\n', results_dir);
end

% S1 significance threshold: parameters with ST below this are candidates
% for fixing to default values during L-BFGS-B optimisation.
ST_threshold = 0.05;

fprintf('  Base sample size N:     %d\n', N);
fprintf('  Patient index:          %d\n', patient_idx);
fprintf('  ST significance threshold: %.2f\n\n', ST_threshold);

%% =========================================================================
%  STEP 1 — GENERATE SOBOL/SALTELLI SAMPLES
% =========================================================================
fprintf('STEP 1: Generating Saltelli sampling matrices...\n');
tic;
[X_all, param_names, param_bounds, n_total, sample_info] = sample_sobol_params(N);
t_sampling = toc;
fprintf('  Sampling complete: %d samples in %.2f s\n\n', n_total, t_sampling);

%% =========================================================================
%  STEP 2 — EVALUATE MODEL FOR ALL SAMPLES
% =========================================================================
fprintf('STEP 2: Evaluating PDA-CoA model for %d parameter combinations...\n', n_total);
fprintf('  This may take 30-90 minutes depending on N.\n');
fprintf('  Progress updates every 50 simulations.\n\n');

tic;
Y = evaluate_model_outputs(X_all, param_names, csv_path, patient_idx);
t_eval = toc;
fprintf('  Total evaluation time: %.1f s (%.1f min)\n\n', t_eval, t_eval/60);

%% =========================================================================
%  STEP 3 — COMPUTE SOBOL INDICES
% =========================================================================
fprintf('STEP 3: Computing Sobol sensitivity indices...\n');

output_names = {'DeltaP_coa_peak', 'DeltaP_coa_mean_sys', 'Q_coa_fraction', 'severity_code'};
n_outputs = length(output_names);

tic;
[S1, ST, S1_ci, ST_ci] = compute_sobol_indices(Y, sample_info);
t_indices = toc;
fprintf('  Index computation complete in %.2f s\n\n', t_indices);

D = sample_info.D;

%% =========================================================================
%  STEP 4 — DISPLAY RESULTS
% =========================================================================
fprintf('STEP 4: Sobol Sensitivity Index Results\n');
fprintf('=========================================================================\n\n');

for q = 1:n_outputs
    fprintf('--- Output: %s ---\n', output_names{q});
    fprintf('  %-18s %8s %12s %8s %12s\n', 'Parameter', 'S1', '95% CI', 'ST', '95% CI');
    fprintf('  %s\n', repmat('-', 1, 62));

    for i = 1:D
        fprintf('  %-18s %8.4f [%5.3f,%5.3f] %8.4f [%5.3f,%5.3f]\n', ...
            param_names{i}, ...
            S1(i, q), S1_ci(i, q, 1), S1_ci(i, q, 2), ...
            ST(i, q), ST_ci(i, q, 1), ST_ci(i, q, 2));
    end
    fprintf('\n');
end

%% =========================================================================
%  STEP 5 — RANK PARAMETERS BY IMPORTANCE
% =========================================================================
fprintf('STEP 5: Parameter Ranking (by mean Total-Order index across outputs)\n');
fprintf('=========================================================================\n');

% Mean ST across all outputs (excluding severity code which is ordinal)
ST_mean_continuous = mean(ST(:, 1:3), 2);   % Mean over DeltaP_peak, DeltaP_mean, Q_frac

[ST_sorted, rank_order] = sort(ST_mean_continuous, 'descend');

fprintf('\n  Rank | %-18s | Mean ST  | Influential?\n', 'Parameter');
fprintf('  %s\n', repmat('-', 1, 55));

influential_params = {};
fixed_params = {};

for r = 1:D
    idx = rank_order(r);
    if ST_sorted(r) >= ST_threshold
        flag = '  ✓  INCLUDE in L-BFGS-B';
        influential_params{end+1} = param_names{idx}; %#ok<SAGROW>
    else
        flag = '  ✗  FIX to default';
        fixed_params{end+1} = param_names{idx}; %#ok<SAGROW>
    end
    fprintf('  %4d | %-18s | %7.4f  |%s\n', r, param_names{idx}, ST_sorted(r), flag);
end

fprintf('\n  Significance threshold: ST ≥ %.2f\n', ST_threshold);
fprintf('  Influential parameters (%d):  %s\n', length(influential_params), strjoin(influential_params, ', '));
fprintf('  Fixed parameters (%d):        %s\n\n', length(fixed_params), strjoin(fixed_params, ', '));

%% =========================================================================
%  STEP 6 — SAVE RESULTS TO CSV
% =========================================================================
fprintf('STEP 6: Saving results to CSV and MAT...\n');

% --- 6a. Full indices table ---
% Build a table with S1, ST, and CIs for each output
rows_cell = {};
for q = 1:n_outputs
    for i = 1:D
        rows_cell{end+1, 1} = output_names{q};          %#ok<SAGROW>
        rows_cell{end, 2}   = param_names{i};
        rows_cell{end, 3}   = S1(i, q);
        rows_cell{end, 4}   = S1_ci(i, q, 1);
        rows_cell{end, 5}   = S1_ci(i, q, 2);
        rows_cell{end, 6}   = ST(i, q);
        rows_cell{end, 7}   = ST_ci(i, q, 1);
        rows_cell{end, 8}   = ST_ci(i, q, 2);
    end
end

T_indices = cell2table(rows_cell, ...
    'VariableNames', {'Output', 'Param
ter', 'S1', 'S1_CI_lower', 'S1_CI_upper', ...
                      'ST', 'ST_CI_lower', 'ST_CI_upper'});

csv_indices_path = fullfile(results_dir, 'sobol_indices.csv');
writetable(T_indices, csv_indices_path);
fprintf('  Saved: %s\n', csv_indices_path);

% --- 6b. Ranking table ---
rank_data = cell(D, 5);
for r = 1:D
    idx = rank_order(r);
    rank_data{r, 1} = r;
    rank_data{r, 2} = param_names{idx};
    rank_data{r, 3} = ST_mean_continuous(idx);
    rank_data{r, 4} = ST_mean_continuous(idx) >= ST_threshold;
    rank_data{r, 5} = S1(idx, 2);  % S1 for DeltaP_mean_sys
end

T_ranking = cell2table(rank_data, ...
    'VariableNames', {'Rank', 'Parameter', 'Mean_ST', 'Influential', 'S1_DeltaP_mean_sys'});

csv_ranking_path = fullfile(results_dir, 'sobol_ranking.csv');
writetable(T_ranking, csv_ranking_path);
fprintf('  Saved: %s\n', csv_ranking_path);

% --- 6c. Save full workspace ---
mat_path = fullfile(results_dir, 'gsa_raw_data.mat');
save(mat_path, 'X_all', 'Y', 'S1', 'ST', 'S1_ci', 'ST_ci', ...
     'param_names', 'param_bounds', 'output_names', 'sample_info', ...
     'N', 'patient_idx', 'influential_params', 'fixed_params', ...
     'ST_threshold', 'rank_order', 'ST_mean_continuous');
fprintf('  Saved: %s\n\n', mat_path);

%% =========================================================================
%  STEP 7 — GENERATE BAR PLOTS
% =========================================================================
fprintf('STEP 7: Generating publication-ready bar plots...\n');

output_labels = {'\DeltaP_{CoA,peak}', '\DeltaP_{CoA,mean-sys}', ...
                 'Q_{CoA}/Q_{total}', 'Severity Code'};

% Colour palette (vibrant, distinguishable)
colors_S1 = [0.22 0.56 0.82];     % Steel blue
colors_ST = [0.89 0.35 0.27];     % Coral red

% --- 7a. First-Order Indices (S1) Bar Plot ---
fig1 = figure('Name', 'Sobol S1 Indices', 'Position', [100 100 1200 700], ...
    'Color', 'w');

for q = 1:n_outputs
    subplot(2, 2, q);

    % Sort by S1 for this output
    [s1_sorted, s_order] = sort(S1(:, q), 'descend');
    ci_lower = S1_ci(s_order, q, 1);
    ci_upper = S1_ci(s_order, q, 2);
    err_neg  = s1_sorted - ci_lower;
    err_pos  = ci_upper - s1_sorted;

    b = barh(1:D, s1_sorted, 'FaceColor', colors_S1, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.85);
    hold on;
    errorbar(s1_sorted, 1:D, err_neg, err_pos, '.', 'horizontal', ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'CapSize', 4);

    % Threshold line
    xline(ST_threshold, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, ...
        'Label', sprintf('Threshold = %.2f', ST_threshold), ...
        'LabelHorizontalAlignment', 'left', 'FontSize', 8);

    set(gca, 'YTick', 1:D, 'YTickLabel', param_names(s_order), ...
        'TickLabelInterpreter', 'none', ...
        'FontSize', 9, 'YDir', 'reverse', 'Box', 'off');
    xlabel('S_1 (First-Order Index)', 'FontSize', 10);
    title(output_labels{q}, 'FontSize', 12, 'FontWeight', 'bold');
    xlim([min(0, min(s1_sorted) - 0.05), max(s1_sorted) + 0.15]);
    grid on;
    set(gca, 'GridAlpha', 0.15);
    hold off;
end

sgtitle('Sobol First-Order Sensitivity Indices (S_1)', ...
    'FontSize', 14, 'FontWeight', 'bold');

fig1_path = fullfile(results_dir, 'sobol_barplot_S1.png');
exportgraphics(fig1, fig1_path, 'Resolution', 300);
fprintf('  Saved: %s\n', fig1_path);

% --- 7b. Total-Order Indices (ST) Bar Plot ---
fig2 = figure('Name', 'Sobol ST Indices', 'Position', [150 100 1200 700], ...
    'Color', 'w');

for q = 1:n_outputs
    subplot(2, 2, q);

    [st_sorted, s_order] = sort(ST(:, q), 'descend');
    ci_lower = ST_ci(s_order, q, 1);
    ci_upper = ST_ci(s_order, q, 2);
    err_neg  = st_sorted - ci_lower;
    err_pos  = ci_upper - st_sorted;

    b = barh(1:D, st_sorted, 'FaceColor', colors_ST, 'EdgeColor', 'none', ...
        'FaceAlpha', 0.85);
    hold on;
    errorbar(st_sorted, 1:D, err_neg, err_pos, '.', 'horizontal', ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'CapSize', 4);

    xline(ST_threshold, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, ...
        'Label', sprintf('Threshold = %.2f', ST_threshold), ...
        'LabelHorizontalAlignment', 'left', 'FontSize', 8);

    set(gca, 'YTick', 1:D, 'YTickLabel', param_names(s_order), ...
        'TickLabelInterpreter', 'none', ...
        'FontSize', 9, 'YDir', 'reverse', 'Box', 'off');
    xlabel('S_T (Total-Order Index)', 'FontSize', 10);
    title(output_labels{q}, 'FontSize', 12, 'FontWeight', 'bold');
    xlim([min(0, min(st_sorted) - 0.05), max(st_sorted) + 0.15]);
    grid on;
    set(gca, 'GridAlpha', 0.15);
    hold off;
end

sgtitle('Sobol Total-Order Sensitivity Indices (S_T)', ...
    'FontSize', 14, 'FontWeight', 'bold');

fig2_path = fullfile(results_dir, 'sobol_barplot_ST.png');
exportgraphics(fig2, fig2_path, 'Resolution', 300);
fprintf('  Saved: %s\n', fig2_path);

% --- 7c. Combined S1 + ST Grouped Bar Plot (for ΔP_CoA_mean_sys) ---
fig3 = figure('Name', 'Sobol S1 vs ST — Mean Systolic Gradient', ...
    'Position', [200 100 900 550], 'Color', 'w');

q_main = 2;  % DeltaP_coa_mean_sys (primary clinical output)

% Use rank_order from Step 5 so plot and printed ranking are consistent
% (both rank by mean ST across the three continuous outputs)
rank_idx = rank_order;
S1_ranked = S1(rank_idx, q_main);
ST_ranked = ST(rank_idx, q_main);

bar_data = [S1_ranked, ST_ranked];
b = barh(1:D, bar_data, 'grouped');
b(1).FaceColor = colors_S1;
b(2).FaceColor = colors_ST;
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';
b(1).FaceAlpha = 0.85;
b(2).FaceAlpha = 0.85;

hold on;
xline(ST_threshold, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, ...
    'Label', sprintf('Threshold = %.2f', ST_threshold), ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 9);

set(gca, 'YTick', 1:D, 'YTickLabel', param_names(rank_idx), ...
    'TickLabelInterpreter', 'none', ...
    'FontSize', 10, 'YDir', 'reverse', 'Box', 'off');
xlabel('Sobol Index', 'FontSize', 12);
title('\DeltaP_{CoA,mean-sys}: First-Order vs Total-Order Sensitivity', ...
    'FontSize', 14, 'FontWeight', 'bold');
legend({'S_1 (First-Order)', 'S_T (Total-Order)'}, 'Location', 'southeast', ...
    'FontSize', 10);
grid on;
set(gca, 'GridAlpha', 0.15);
hold off;

fig3_path = fullfile(results_dir, 'sobol_barplot_combined.png');
exportgraphics(fig3, fig3_path, 'Resolution', 300);
fprintf('  Saved: %s\n\n', fig3_path);

%% =========================================================================
%  STEP 8 — SUMMARY AND RECOMMENDATIONS
% =========================================================================
fprintf('=========================================================================\n');
fprintf('   GSA SUMMARY AND RECOMMENDATIONS FOR L-BFGS-B OPTIMISATION\n');
fprintf('=========================================================================\n\n');

fprintf('  ┌────────────────────────────────────────────────────────────────┐\n');
fprintf('  │  INFLUENTIAL PARAMETERS (include in L-BFGS-B)                 │\n');
fprintf('  │  ST ≥ %.2f across continuous clinical outputs                 │\n', ST_threshold);
fprintf('  ├────────────────────────────────────────────────────────────────┤\n');
for j = 1:length(influential_params)
    idx_p = find(strcmp(param_names, influential_params{j}));
    fprintf('  │  %d. %-18s  ST_mean = %.4f                      │\n', ...
        j, influential_params{j}, ST_mean_continuous(idx_p));
end
fprintf('  └────────────────────────────────────────────────────────────────┘\n\n');

fprintf('  ┌────────────────────────────────────────────────────────────────┐\n');
fprintf('  │  FIXED PARAMETERS (set to default/literature values)          │\n');
fprintf('  │  ST < %.2f — negligible influence on clinical outputs         │\n', ST_threshold);
fprintf('  ├────────────────────────────────────────────────────────────────┤\n');
if isempty(fixed_params)
    fprintf('  │  (none — all parameters are influential)                     │\n');
else
    for j = 1:length(fixed_params)
        idx_p = find(strcmp(param_names, fixed_params{j}));
        fprintf('  │  %d. %-18s  ST_mean = %.4f  → FIX to default   │\n', ...
            j, fixed_params{j}, ST_mean_continuous(idx_p));
    end
end
fprintf('  └────────────────────────────────────────────────────────────────┘\n\n');

fprintf('  Total computation time: %.1f min\n', (t_sampling + t_eval + t_indices) / 60);
fprintf('  Results saved to: %s/\n', results_dir);
fprintf('  Figures saved as PNG (300 dpi).\n');
fprintf('\n  GSA complete. Ready for L-BFGS-B parameter optimisation.\n');
fprintf('=========================================================================\n');
