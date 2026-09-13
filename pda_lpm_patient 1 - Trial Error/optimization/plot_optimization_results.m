function plot_optimization_results(baseline_outputs, opt_outputs, clinical, ...
                                    x_history, J_history, param_names, results_dir)
% PLOT_OPTIMIZATION_RESULTS
% -----------------------------------------------------------------------
% Generates and saves publication-ready plots for the PDA-only Patient 1
% optimization (MAP + SV objective only).
%
% Two figures are produced and saved to results_dir:
%
%   Figure 1 — Objective function convergence
%     Plots J (weighted error) vs optimizer iteration number.
%
%   Figure 2 — Primary Optimization Outputs: Before vs After
%     Grouped bar chart comparing MAP and SV between clinical
%     measurements and both simulations (pre- and post-optimization).
%     SBP, DBP, dP_PDA, and dP_CoA are NOT shown here because they
%     are not direct optimization targets for this configuration.
%
% No CoA clinical target plot is produced.  dP_CoA = 0 in the CSV
% reflects a missing measurement, not a clinical target of 0 mmHg.
%
% INPUTS:
%   baseline_outputs - struct from objective function at x0 (pre-opt)
%   opt_outputs      - struct from objective function at x_opt (post-opt)
%   clinical         - clinical measurement struct
%   x_history        - (n_iter × D_opt) matrix of parameter values per iter
%   J_history        - (n_iter × 1) vector of J values per iteration
%   param_names      - cell array of optimized parameter names
%   results_dir      - directory to save PNG files
%
% AUTHOR:   Optimization Extension — Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  2.0  — restricted to MAP+SV for PDA-only patient
% -----------------------------------------------------------------------

% Ensure output directory exists
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% -------------------------------------------------------------------------
% Colour palette
% -------------------------------------------------------------------------
C_clinical  = [0.20, 0.53, 0.74];   % Steel blue — clinical reference
C_baseline  = [0.85, 0.50, 0.20];   % Burnt orange — pre-optimization
C_optimized = [0.22, 0.66, 0.40];   % Emerald green — post-optimization
C_objective = [0.60, 0.20, 0.60];   % Purple — objective trace

% =========================================================================
%  FIGURE 1 — Objective Function Convergence
% =========================================================================
fig1 = figure('Name', 'Objective Convergence', ...
    'Position', [80 80 850 450], 'Color', 'w');

% Remove any failed (NaN/Inf) entries from J_history for clean plotting
valid_J  = isfinite(J_history);
iter_vec = find(valid_J);
J_valid  = J_history(valid_J);

subplot(1, 2, 1);
semilogy(iter_vec, J_valid, '-o', ...
    'Color', C_objective, 'LineWidth', 1.8, 'MarkerSize', 4, ...
    'MarkerFaceColor', C_objective);
xlabel('Iteration', 'FontSize', 11);
ylabel('Objective J (log scale)', 'FontSize', 11);
title('Convergence: J vs Iteration', 'FontSize', 13, 'FontWeight', 'bold');
grid on; set(gca, 'GridAlpha', 0.18, 'Box', 'off', 'FontSize', 10);
hold on;
% Mark initial and final values
if ~isempty(J_valid)
    plot(iter_vec(1),   J_valid(1),   'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', C_baseline,  'MarkerEdgeColor', 'k', 'LineWidth', 1);
    plot(iter_vec(end), J_valid(end), 's', 'MarkerSize', 10, ...
        'MarkerFaceColor', C_optimized, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    legend({'J trace', 'Initial J', 'Final J'}, 'Location', 'northeast', 'FontSize', 9);
end
hold off;

% Parameter trajectories (normalized to initial value)
subplot(1, 2, 2);
hold on;
D_opt = size(x_history, 2);
cmap  = lines(D_opt);
for k = 1:D_opt
    x_k = x_history(valid_J, k);
    if x_history(1, k) ~= 0
        x_k_norm = x_k / x_history(1, k);   % Normalize to initial value
    else
        x_k_norm = x_k;
    end
    plot(iter_vec, x_k_norm, '-', 'Color', cmap(k, :), 'LineWidth', 1.5);
end
xlabel('Iteration', 'FontSize', 11);
ylabel('Parameter value / Initial value', 'FontSize', 11);
title('Parameter Evolution (Normalized)', 'FontSize', 13, 'FontWeight', 'bold');
legend(param_names, 'Location', 'best', 'FontSize', 8, 'Interpreter', 'none');
yline(1, '--k', 'LineWidth', 0.8, 'Alpha', 0.4);
grid on; set(gca, 'GridAlpha', 0.18, 'Box', 'off', 'FontSize', 10);
hold off;

sgtitle('L-BFGS-B Optimization Convergence', 'FontSize', 14, 'FontWeight', 'bold');

fig1_path = fullfile(results_dir, 'opt_convergence.png');
exportgraphics(fig1, fig1_path, 'Resolution', 300);
fprintf('  Saved: %s\n', fig1_path);

% =========================================================================
%  FIGURE 2 — Primary Optimization Outputs: Before vs After
% =========================================================================
% Only MAP and SV are shown — these are the direct clinical targets for
% the PDA-only Patient 1 configuration.  SBP, DBP, dP_PDA, and dP_CoA
% are not direct targets in this 3-parameter setup and are not plotted.
fig2 = figure('Name', 'Primary Optimization Outputs', ...
    'Position', [100 100 700 500], 'Color', 'w');

target_labels = {'MAP (mmHg)', 'SV (mL)'};
clinical_vals = [clinical.P_ao_mean_mmHg; clinical.SV_mL];

baseline_vals = [safe_get(baseline_outputs, 'P_ao_mean'); ...
                 safe_get(baseline_outputs, 'SV_lv')];
opt_vals      = [safe_get(opt_outputs, 'P_ao_mean'); ...
                 safe_get(opt_outputs, 'SV_lv')];

n_targets = length(target_labels);
bar_data  = [clinical_vals, baseline_vals, opt_vals];

b = bar(1:n_targets, bar_data, 'grouped');
b(1).FaceColor = C_clinical;   b(1).EdgeColor = 'none'; b(1).FaceAlpha = 0.9;
b(2).FaceColor = C_baseline;   b(2).EdgeColor = 'none'; b(2).FaceAlpha = 0.8;
b(3).FaceColor = C_optimized;  b(3).EdgeColor = 'none'; b(3).FaceAlpha = 0.9;

set(gca, 'XTick', 1:n_targets, 'XTickLabel', target_labels, ...
    'FontSize', 11, 'Box', 'off', 'GridAlpha', 0.15);
ylabel('Value', 'FontSize', 12);
title('Primary Optimization Outputs: Before vs After Optimization', ...
    'FontSize', 13, 'FontWeight', 'bold');
legend({'Clinical Reference', 'Baseline (Pre-Opt)', 'Optimized (Post-Opt)'}, ...
    'Location', 'northwest', 'FontSize', 10);
grid on;

% Annotate % error (normalized to max(|clin|, floor) — same as objective)
floors = [1.0; 0.1];
for i = 1:n_targets
    norm_denom = max(abs(clinical_vals(i)), floors(i));
    if isfinite(baseline_vals(i)) && isfinite(opt_vals(i))
        err_base = 100 * abs(baseline_vals(i) - clinical_vals(i)) / norm_denom;
        err_opt  = 100 * abs(opt_vals(i)      - clinical_vals(i)) / norm_denom;
        text(i + 0.22, opt_vals(i) * 1.02, sprintf('%.1f%%', err_opt), ...
            'FontSize', 8, 'Color', C_optimized, 'HorizontalAlignment', 'center');
        text(i - 0.02, baseline_vals(i) * 1.02, sprintf('%.1f%%', err_base), ...
            'FontSize', 8, 'Color', C_baseline,  'HorizontalAlignment', 'center');
    end
end

fig2_path = fullfile(results_dir, 'opt_clinical_comparison.png');
exportgraphics(fig2, fig2_path, 'Resolution', 300);
fprintf('  Saved: %s\n', fig2_path);
% Note: CoA summary figure (opt_coa_summary.png) is not produced for
% Patient 1 because CoA gradient is not a measured clinical target.

end


%% -------------------------------------------------------------------------
%  LOCAL HELPER: safe_get
%  Returns field value from struct, or NaN if field absent
% -------------------------------------------------------------------------
function val = safe_get(s, fname)
    if isstruct(s) && isfield(s, fname)
        val = s.(fname);
        if ~isnumeric(val) || isempty(val)
            val = NaN;
        end
    else
        val = NaN;
    end
end

