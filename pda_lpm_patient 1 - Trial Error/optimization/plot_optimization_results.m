function plot_optimization_results(baseline_outputs, opt_outputs, clinical, ...
    x_history, J_history, param_names, results_dir)
% PLOT_OPTIMIZATION_RESULTS
% -------------------------------------------------------------------------
% Creates publication-ready plots for PDA patient-specific optimization.
%
% Figure 1: objective convergence and parameter trajectories.
% Figure 2: clinical, baseline, and optimized values for every available
%           calibration target.
% Figure 3: absolute percentage error before and after optimization.
%
% Core targets: MAP, SV, SBP, DBP, and peak PDA pressure gradient.
% Optional targets (shown only when measured): CoA gradient, Qp/Qs,
% LV ejection fraction, and LV end-diastolic volume.
% -------------------------------------------------------------------------

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% Colour palette
C_clinical  = [0.20, 0.53, 0.74];
C_baseline  = [0.85, 0.50, 0.20];
C_optimized = [0.22, 0.66, 0.40];
C_objective = [0.60, 0.20, 0.60];

%% ========================================================================
% FIGURE 1: OBJECTIVE AND PARAMETER CONVERGENCE
% =========================================================================
n_hist = min(size(x_history, 1), numel(J_history));
if n_hist < 1
    warning('PLOT_OPTIMIZATION_RESULTS:NoHistory', ...
        'No optimization history was available for plotting.');
else
    x_hist = x_history(1:n_hist, :);
    J_hist = J_history(1:n_hist);
    iter_all = (0:n_hist-1)';
    valid_J = isfinite(J_hist);
    iter_vec = iter_all(valid_J);
    J_valid = J_hist(valid_J);
    x_valid = x_hist(valid_J, :);

    fig1 = figure('Name', 'Optimization Convergence', ...
        'Position', [80 80 1100 470], 'Color', 'w');

    subplot(1, 2, 1);
    semilogy(iter_vec, max(J_valid, eps), '-o', ...
        'Color', C_objective, 'LineWidth', 1.8, 'MarkerSize', 4, ...
        'MarkerFaceColor', C_objective);
    xlabel('Iteration', 'FontSize', 11);
    ylabel('Objective J (log scale)', 'FontSize', 11);
    title('Objective Convergence', 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
    set(gca, 'GridAlpha', 0.18, 'Box', 'off', 'FontSize', 10);
    hold on;
    if ~isempty(J_valid)
        plot(iter_vec(1), max(J_valid(1), eps), 'o', 'MarkerSize', 9, ...
            'MarkerFaceColor', C_baseline, 'MarkerEdgeColor', 'k');
        plot(iter_vec(end), max(J_valid(end), eps), 's', 'MarkerSize', 9, ...
            'MarkerFaceColor', C_optimized, 'MarkerEdgeColor', 'k');
        legend({'J trace', 'Initial J', 'Final J'}, ...
            'Location', 'northeast', 'FontSize', 9);
    end
    hold off;

    subplot(1, 2, 2);
    hold on;
    D_opt = min(size(x_valid, 2), numel(param_names));
    cmap = lines(max(D_opt, 1));
    for k = 1:D_opt
        x_k = x_valid(:, k);
        x_initial = x_hist(1, k);
        if isfinite(x_initial) && abs(x_initial) > eps
            x_k = x_k ./ x_initial;
        end
        plot(iter_vec, x_k, '-', 'Color', cmap(k, :), 'LineWidth', 1.5);
    end
    yline(1, '--k', 'LineWidth', 0.8, 'Alpha', 0.4);
    xlabel('Iteration', 'FontSize', 11);
    ylabel('Parameter / Initial Value', 'FontSize', 11);
    title('Parameter Evolution', 'FontSize', 13, 'FontWeight', 'bold');
    if D_opt > 0
        legend(param_names(1:D_opt), 'Location', 'best', ...
            'FontSize', 8, 'Interpreter', 'none');
    end
    grid on;
    set(gca, 'GridAlpha', 0.18, 'Box', 'off', 'FontSize', 10);
    hold off;

    sgtitle('Bounded Optimization Convergence', ...
        'FontSize', 14, 'FontWeight', 'bold');

    fig1_path = fullfile(results_dir, 'opt_convergence.png');
    exportgraphics(fig1, fig1_path, 'Resolution', 300);
    fprintf('  Saved: %s\n', fig1_path);
end

%% ========================================================================
% BUILD TARGET LIST
% =========================================================================
target_labels = {'MAP', 'SV', 'SBP', 'DBP', 'PDA Gradient'};
target_units  = {'mmHg', 'mL', 'mmHg', 'mmHg', 'mmHg'};

clinical_vals = [ ...
    safe_get(clinical, 'P_ao_mean_mmHg'); ...
    safe_get(clinical, 'SV_mL'); ...
    safe_get(clinical, 'P_ao_sys_mmHg'); ...
    safe_get(clinical, 'P_ao_dia_mmHg'); ...
    safe_get(clinical, 'dP_pda_mmHg')];

baseline_vals = [ ...
    safe_get(baseline_outputs, 'P_ao_mean'); ...
    safe_get(baseline_outputs, 'SV_lv'); ...
    safe_get(baseline_outputs, 'P_ao_sys'); ...
    safe_get(baseline_outputs, 'P_ao_dia'); ...
    safe_get(baseline_outputs, 'dP_PDA_peak')];

opt_vals = [ ...
    safe_get(opt_outputs, 'P_ao_mean'); ...
    safe_get(opt_outputs, 'SV_lv'); ...
    safe_get(opt_outputs, 'P_ao_sys'); ...
    safe_get(opt_outputs, 'P_ao_dia'); ...
    safe_get(opt_outputs, 'dP_PDA_peak')];

normalization_floors = [1.0; 0.1; 1.0; 1.0; 0.5];

% Optional measured CoA gradient
clin_coa = safe_get(clinical, 'dP_coa_mmHg');
if isfinite(clin_coa) && clin_coa > 0
    target_labels{end+1} = 'CoA Gradient'; %#ok<AGROW>
    target_units{end+1} = 'mmHg'; %#ok<AGROW>
    clinical_vals(end+1,1) = clin_coa; %#ok<AGROW>
    baseline_vals(end+1,1) = safe_get(baseline_outputs, 'DeltaP_coa_peak'); %#ok<AGROW>
    opt_vals(end+1,1) = safe_get(opt_outputs, 'DeltaP_coa_peak'); %#ok<AGROW>
    normalization_floors(end+1,1) = 1.0; %#ok<AGROW>
end

% Optional measured Qp/Qs
clin_qpqs = safe_get(clinical, 'Qp_Qs');
if isfinite(clin_qpqs) && clin_qpqs > 0
    target_labels{end+1} = 'Qp/Qs'; %#ok<AGROW>
    target_units{end+1} = 'ratio'; %#ok<AGROW>
    clinical_vals(end+1,1) = clin_qpqs; %#ok<AGROW>
    baseline_vals(end+1,1) = safe_get(baseline_outputs, 'Qp_Qs'); %#ok<AGROW>
    opt_vals(end+1,1) = safe_get(opt_outputs, 'Qp_Qs'); %#ok<AGROW>
    normalization_floors(end+1,1) = 0.1; %#ok<AGROW>
end

% Optional measured LV ejection fraction
clin_ef = safe_get(clinical, 'EF_lv_pct');
base_ef = safe_get(baseline_outputs, 'EF_lv');
opt_ef  = safe_get(opt_outputs, 'EF_lv');
if isfinite(clin_ef) && clin_ef > 0
    if clin_ef > 1
        if isfinite(base_ef) && base_ef <= 1, base_ef = 100 * base_ef; end
        if isfinite(opt_ef)  && opt_ef  <= 1, opt_ef  = 100 * opt_ef;  end
    end
    target_labels{end+1} = 'LV EF'; %#ok<AGROW>
    target_units{end+1} = '%'; %#ok<AGROW>
    clinical_vals(end+1,1) = clin_ef; %#ok<AGROW>
    baseline_vals(end+1,1) = base_ef; %#ok<AGROW>
    opt_vals(end+1,1) = opt_ef; %#ok<AGROW>
    normalization_floors(end+1,1) = 1.0; %#ok<AGROW>
end

% Optional measured LV end-diastolic volume
clin_edv = safe_get(clinical, 'EDV_lv_mL');
if isfinite(clin_edv) && clin_edv > 0
    target_labels{end+1} = 'LV EDV'; %#ok<AGROW>
    target_units{end+1} = 'mL'; %#ok<AGROW>
    clinical_vals(end+1,1) = clin_edv; %#ok<AGROW>
    baseline_vals(end+1,1) = safe_get(baseline_outputs, 'EDV_lv'); %#ok<AGROW>
    opt_vals(end+1,1) = safe_get(opt_outputs, 'EDV_lv'); %#ok<AGROW>
    normalization_floors(end+1,1) = 0.1; %#ok<AGROW>
end

% Do not plot targets whose clinical reference is unavailable.
valid_targets = isfinite(clinical_vals);
target_labels = target_labels(valid_targets);
target_units = target_units(valid_targets);
clinical_vals = clinical_vals(valid_targets);
baseline_vals = baseline_vals(valid_targets);
opt_vals = opt_vals(valid_targets);
normalization_floors = normalization_floors(valid_targets);
n_targets = numel(target_labels);

%% ========================================================================
% FIGURE 2: CLINICAL TARGETS BEFORE AND AFTER OPTIMIZATION
% =========================================================================
if n_targets > 0
    n_cols = min(3, n_targets);
    n_rows = ceil(n_targets / n_cols);
    fig2_height = max(440, 320 * n_rows);
    fig2 = figure('Name', 'Clinical Target Comparison', ...
        'Position', [100 60 1150 fig2_height], 'Color', 'w');
    tl = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    for i = 1:n_targets
        ax = nexttile;
        vals = [clinical_vals(i), baseline_vals(i), opt_vals(i)];
        b = bar(ax, 1:3, vals, 0.70, 'FaceColor', 'flat', ...
            'EdgeColor', 'none');
        b.CData = [C_clinical; C_baseline; C_optimized];
        set(ax, 'XTick', 1:3, 'XTickLabel', ...
            {'Clinical', 'Pre-Opt', 'Post-Opt'}, ...
            'FontSize', 9, 'Box', 'off', 'GridAlpha', 0.15);
        ylabel(ax, target_units{i}, 'FontSize', 10);
        title(ax, target_labels{i}, 'FontSize', 11, 'FontWeight', 'bold');
        grid(ax, 'on');

        finite_vals = vals(isfinite(vals));
        if ~isempty(finite_vals)
            ymax = max(finite_vals);
            ymin = min([0, finite_vals]);
            span = max(ymax - ymin, max(abs(finite_vals)) * 0.2);
            if span <= 0, span = 1; end
            ylim(ax, [ymin - 0.05*span, ymax + 0.28*span]);
        end

        denom = max(abs(clinical_vals(i)), normalization_floors(i));
        if isfinite(baseline_vals(i))
            err_pre = 100 * abs(baseline_vals(i) - clinical_vals(i)) / denom;
            text(ax, 2, baseline_vals(i), sprintf('  %.1f%% ', err_pre), ...
                'Color', C_baseline, 'FontSize', 8, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
        if isfinite(opt_vals(i))
            err_post = 100 * abs(opt_vals(i) - clinical_vals(i)) / denom;
            text(ax, 3, opt_vals(i), sprintf('  %.1f%% ', err_post), ...
                'Color', C_optimized, 'FontSize', 8, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end

    title(tl, 'Clinical Targets: Before and After Optimization', ...
        'FontSize', 14, 'FontWeight', 'bold');

    fig2_path = fullfile(results_dir, 'opt_clinical_comparison.png');
    exportgraphics(fig2, fig2_path, 'Resolution', 300);
    fprintf('  Saved: %s\n', fig2_path);
end

end

%% ------------------------------------------------------------------------
% LOCAL HELPER
% -------------------------------------------------------------------------
function val = safe_get(s, field_name)
if isstruct(s) && isfield(s, field_name)
    val = s.(field_name);
    if ~isnumeric(val) || ~isscalar(val) || isempty(val)
        val = NaN;
    end
else
    val = NaN;
end
end