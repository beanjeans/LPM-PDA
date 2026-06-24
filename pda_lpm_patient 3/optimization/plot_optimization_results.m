function plot_optimization_results(baseline_outputs, opt_outputs, clinical, ...
                                    x_history, J_history, param_names, results_dir)
% PLOT_OPTIMIZATION_RESULTS
% -----------------------------------------------------------------------
% Generates and saves publication-ready plots comparing the cardiovascular
% model simulation before and after L-BFGS-B parameter optimization.
%
% Three figures are produced and saved to results_dir:
%
%   Figure 1 — Objective function convergence
%     Plots J (weighted error) vs optimizer iteration number, with a
%     log-scale y-axis for clarity.
%
%   Figure 2 — Clinical vs Simulated: Before vs After
%     Grouped bar chart comparing key haemodynamic targets (MAP, SBP, DBP,
%     SV, PDA gradient) between clinical measurements and both simulations.
%
%   Figure 3 — CoA Clinical Output Summary
%     Displays the four key CoA outputs (ΔP_peak, ΔP_mean_sys,
%     Q_coa/Q_total, predicted severity) before and after optimization.
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
% VERSION:  1.0
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
%  Only shows MAP, SV, and dP_CoA_peak — the direct clinical targets.
% =========================================================================
fig2 = figure('Name', 'Primary Optimization Outputs', ...
    'Position', [100 100 800 520], 'Color', 'w');

% Direct optimization targets only: MAP, SV, dP_CoA_peak
target_labels  = {'MAP (mmHg)', 'SV (mL)', 'dP_{CoA,peak} (mmHg)'};
clinical_vals  = [
    clinical.P_ao_mean_mmHg
    clinical.SV_mL
    clinical.dP_coa_mmHg   % measured CoA target (dPCoA = 4.9 mmHg for Patient 3)
];

% Safe extraction with fallback to NaN if field missing
baseline_vals = [
    safe_get(baseline_outputs, 'P_ao_mean')
    safe_get(baseline_outputs, 'SV_lv')
    safe_get(baseline_outputs, 'DeltaP_coa_peak')
];

opt_vals = [
    safe_get(opt_outputs, 'P_ao_mean')
    safe_get(opt_outputs, 'SV_lv')
    safe_get(opt_outputs, 'DeltaP_coa_peak')
];

n_targets = length(target_labels);
bar_data  = [clinical_vals, baseline_vals, opt_vals];

b = bar(1:n_targets, bar_data, 'grouped');
b(1).FaceColor = C_clinical;   b(1).EdgeColor = 'none'; b(1).FaceAlpha = 0.9;
b(2).FaceColor = C_baseline;   b(2).EdgeColor = 'none'; b(2).FaceAlpha = 0.8;
b(3).FaceColor = C_optimized;  b(3).EdgeColor = 'none'; b(3).FaceAlpha = 0.9;

set(gca, 'XTick', 1:n_targets, 'XTickLabel', target_labels, ...
    'FontSize', 10, 'Box', 'off', 'GridAlpha', 0.15);
ylabel('Value', 'FontSize', 12);
title('Primary Optimization Outputs: Before vs After Optimization', ...
    'FontSize', 13, 'FontWeight', 'bold');
legend({'Clinical Reference', 'Baseline (Pre-Opt)', 'Optimized (Post-Opt)'}, ...
    'Location', 'northwest', 'FontSize', 10);
grid on;
xtickangle(20);

% Annotate % error for MAP, SV, and dP_CoA_peak only
for i = 1:n_targets
    if isfinite(clinical_vals(i)) && clinical_vals(i) > 0
        err_base = 100 * abs(baseline_vals(i) - clinical_vals(i)) / clinical_vals(i);
        err_opt  = 100 * abs(opt_vals(i)      - clinical_vals(i)) / clinical_vals(i);
        if isfinite(opt_vals(i))
            text(i + 0.22, opt_vals(i) * 1.02, sprintf('%.1f%%', err_opt), ...
                'FontSize', 7, 'Color', C_optimized, 'HorizontalAlignment', 'center');
        end
        if isfinite(baseline_vals(i))
            text(i - 0.02, baseline_vals(i) * 1.02, sprintf('%.1f%%', err_base), ...
                'FontSize', 7, 'Color', C_baseline, 'HorizontalAlignment', 'center');
        end
    end
end

% Annotate mild-zone threshold on the dP_CoA_peak bar (index 3)
hold on;
coa_mild_limit = 10.0;  % mmHg — mild-zone upper limit
x_coa = 3;              % dP_CoA_peak is the 3rd bar group
plot([x_coa - 0.4, x_coa + 0.4], [coa_mild_limit, coa_mild_limit], ...
    '--', 'Color', [0.55 0.20 0.55], 'LineWidth', 1.5);
text(x_coa + 0.42, coa_mild_limit, 'mild zone \leq 10 mmHg', ...
    'FontSize', 7, 'Color', [0.55 0.20 0.55], ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
text(x_coa, -max(clinical_vals) * 0.08, 'dP_{CoA} objective: zero penalty if \leq 10 mmHg', ...
    'FontSize', 7, 'Color', [0.55 0.20 0.55], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Units', 'data');
hold off;

fig2_path = fullfile(results_dir, 'opt_clinical_comparison.png');
exportgraphics(fig2, fig2_path, 'Resolution', 300);
fprintf('  Saved: %s\n', fig2_path);

% =========================================================================
%  FIGURE 3 — CoA Clinical Output Summary
% =========================================================================
fig3 = figure('Name', 'CoA Output Summary', ...
    'Position', [120 80 1000 420], 'Color', 'w');

% --- Subplot 1: Pressure gradient comparison ---
subplot(1, 3, 1);
coa_labels = {'\DeltaP_{peak}', '\DeltaP_{mean-sys}'};
coa_base   = [safe_get(baseline_outputs, 'DeltaP_coa_peak'), ...
              safe_get(baseline_outputs, 'DeltaP_coa_mean_sys')];
coa_opt    = [safe_get(opt_outputs, 'DeltaP_coa_peak'), ...
              safe_get(opt_outputs, 'DeltaP_coa_mean_sys')];

if isfield(clinical, 'dP_coa_mmHg') && clinical.dP_coa_mmHg > 0
    coa_clin = [clinical.dP_coa_mmHg, clinical.dP_coa_mmHg * 0.7];  % Approx mean-sys
else
    coa_clin = [NaN, NaN];
end

bar_coa = [coa_clin', coa_base', coa_opt'];
b3 = bar(1:2, bar_coa, 'grouped');
b3(1).FaceColor = C_clinical;  b3(1).EdgeColor = 'none'; b3(1).FaceAlpha = 0.9;
b3(2).FaceColor = C_baseline;  b3(2).EdgeColor = 'none'; b3(2).FaceAlpha = 0.8;
b3(3).FaceColor = C_optimized; b3(3).EdgeColor = 'none'; b3(3).FaceAlpha = 0.9;
set(gca, 'XTick', 1:2, 'XTickLabel', coa_labels, 'FontSize', 10, 'Box', 'off');
ylabel('Pressure Gradient (mmHg)', 'FontSize', 11);
title('\DeltaP_{CoA}', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Clinical (Echo)', 'Pre-Opt', 'Post-Opt'}, 'Location', 'northwest', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.15);

% Reference lines for severity thresholds
yline(20, '--', 'Color', [0.7 0.7 0.1], 'LineWidth', 1.2, ...
    'Label', 'Mild/Mod (20 mmHg)', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);
yline(40, '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.2, ...
    'Label', 'Mod/Sev (40 mmHg)', 'LabelHorizontalAlignment', 'right', 'FontSize', 8);

% --- Subplot 2: Q_coa fraction ---
subplot(1, 3, 2);
q_base = safe_get(baseline_outputs, 'Q_coa_fraction');
q_opt  = safe_get(opt_outputs,      'Q_coa_fraction');
bq = bar([1, 2], [q_base, q_opt], 0.5, 'FaceColor', 'flat');
bq.CData(1, :) = C_baseline;
bq.CData(2, :) = C_optimized;
bq.EdgeColor = 'none';
bq.FaceAlpha = 0.85;
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Pre-Opt', 'Post-Opt'}, ...
    'FontSize', 10, 'Box', 'off');
ylabel('Q_{CoA} / Q_{total}', 'FontSize', 11);
title('CoA Flow Fraction', 'FontSize', 13, 'FontWeight', 'bold');
ylim([0, 1]);
yline(0.5, '--k', 'LineWidth', 1.0, 'Alpha', 0.3);
grid on; set(gca, 'GridAlpha', 0.15);
text(1, q_base + 0.03, sprintf('%.2f', q_base), 'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', C_baseline);
text(2, q_opt  + 0.03, sprintf('%.2f', q_opt),  'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', C_optimized);

% --- Subplot 3: Severity text summary ---
subplot(1, 3, 3);
axis off;
sev_base = severity_to_num(safe_get_str(baseline_outputs, 'predicted_CoA_severity'));
sev_opt  = severity_to_num(safe_get_str(opt_outputs,      'predicted_CoA_severity'));

% Colour-coded severity labels
sev_colors = {[0.20 0.60 0.20], [0.85 0.70 0.10], [0.80 0.15 0.15]};
sev_strs   = {'MILD', 'MODERATE', 'SEVERE'};

text(0.5, 0.88, 'Predicted CoA Severity', 'HorizontalAlignment', 'center', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Units', 'normalized');

text(0.5, 0.72, 'Pre-Optimization:', 'HorizontalAlignment', 'center', ...
    'FontSize', 11, 'Units', 'normalized', 'Color', [0.4 0.4 0.4]);
if sev_base > 0
    text(0.5, 0.56, sev_strs{sev_base}, 'HorizontalAlignment', 'center', ...
        'FontSize', 18, 'FontWeight', 'bold', 'Units', 'normalized', ...
        'Color', sev_colors{sev_base});
else
    text(0.5, 0.56, safe_get_str(baseline_outputs, 'predicted_CoA_severity'), ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'Units', 'normalized');
end

text(0.5, 0.40, 'Post-Optimization:', 'HorizontalAlignment', 'center', ...
    'FontSize', 11, 'Units', 'normalized', 'Color', [0.4 0.4 0.4]);
if sev_opt > 0
    text(0.5, 0.24, sev_strs{sev_opt}, 'HorizontalAlignment', 'center', ...
        'FontSize', 18, 'FontWeight', 'bold', 'Units', 'normalized', ...
        'Color', sev_colors{sev_opt});
else
    text(0.5, 0.24, safe_get_str(opt_outputs, 'predicted_CoA_severity'), ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'Units', 'normalized');
end

sgtitle('CoA Clinical Output: Pre vs Post Optimization', ...
    'FontSize', 14, 'FontWeight', 'bold');

fig3_path = fullfile(results_dir, 'opt_coa_summary.png');
exportgraphics(fig3, fig3_path, 'Resolution', 300);
fprintf('  Saved: %s\n', fig3_path);

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

%% -------------------------------------------------------------------------
%  LOCAL HELPER: safe_get_str
%  Returns string field from struct, or empty string if absent
% -------------------------------------------------------------------------
function val = safe_get_str(s, fname)
    if isstruct(s) && isfield(s, fname)
        val = s.(fname);
        if ~ischar(val) && ~isstring(val)
            val = 'N/A';
        end
    else
        val = 'N/A';
    end
end

%% -------------------------------------------------------------------------
%  LOCAL HELPER: severity_to_num
%  Converts severity string to index 1/2/3 for plotting
% -------------------------------------------------------------------------
function n = severity_to_num(sev_str)
    switch lower(char(sev_str))
        case 'mild',     n = 1;
        case 'moderate', n = 2;
        case 'severe',   n = 3;
        otherwise,       n = 0;
    end
end
