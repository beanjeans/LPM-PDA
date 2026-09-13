function plot_results(t_sol_pda, X_sol_pda, results_pda, ...
                      coa_scenarios, patient_id, params_pda)
% PLOT_RESULTS
% -----------------------------------------------------------------------
% Generates publication-ready figures for a single patient.
%
% PDA-only patients (P01, P02 — coa_scenarios empty):
%   Fig 1: Pressure waveforms (LV, Ao, PA) — PDA baseline
%   Fig 2: LV P-V loop — PDA baseline
%   Fig 3: Primary outputs comparison: MAP and SV (clinical vs simulated)
%          Title: "Primary LPM Outputs: MAP and SV"
%
% PDA-CoA patients (P03 — coa_scenarios non-empty):
%   Fig 1: Pressure waveforms (LV, Ao, PA) — PDA baseline
%   Fig 2: LV P-V loop — PDA baseline
%   Fig 3: Trans-CoA gradient traces — all CoA scenarios
%   Fig 4: Bar chart comparison: MAP, CO, SW_lv, DeltaP_coa across scenarios
%
% INPUTS:
%   t_sol_pda      - time vector for PDA-only solution              [s]
%   X_sol_pda      - state matrix for PDA-only solution
%   results_pda    - indices struct from compute_clinical_indices (PDA)
%   coa_scenarios  - cell array of CoA scenario structs (may be empty)
%   patient_id     - patient identifier string                      [-]
%   params_pda     - PDA parameter struct (required; used for state indices
%                    and timing — avoids crash when coa_scenarios is empty)
%
% REFERENCES:
%   [1] Guardrails §11 — publication-ready plot standards
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  2.0  — disease-mode routing; PDA-only primary outputs figure
% -----------------------------------------------------------------------

% Font settings (Guardrail §11.1)
font_name  = 'Arial';
font_sz_ax = 10;
font_sz_lb = 11;
font_sz_ti = 12;

% Colour palette (colorblind-friendly)
col_lv   = [0.00, 0.45, 0.70];   % Blue  — LV
col_ao   = [0.80, 0.40, 0.00];   % Orange — Aorta
col_pa   = [0.00, 0.62, 0.45];   % Green — PA
col_pda  = [0.94, 0.89, 0.26];   % Yellow — PDA shunt
col_coa  = [0.84, 0.37, 0.00];   % Red-orange — CoA gradient

coa_line_styles = {'-', '--', '-.'};

idx_pda  = results_pda.model;   % shorthand (model outputs, not state indices)

%% =========================================================================
%  FIGURE 1: PDA-BASELINE PRESSURE WAVEFORMS
% =========================================================================
% params_pda is passed as 6th argument; fall back to coa_scenarios if omitted
if nargin < 6 || isempty(params_pda)
    if ~isempty(coa_scenarios)
        params_pda = coa_scenarios{1}.params;
    else
        error('plot_results: params_pda (6th argument) is required when coa_scenarios is empty.');
    end
end
pda_idx = params_pda.idx;

fig1 = figure('Units', 'centimeters', 'Position', [0 0 22 14], 'Color', 'w');
sgtitle(sprintf('PDA Baseline Pressure Waveforms — Patient %s', patient_id), ...
    'FontSize', font_sz_ti, 'FontName', font_name, 'FontWeight', 'bold');

% Extract last cycle
T        = params_pda.T_cardiac;
mask_cyc = t_sol_pda >= (t_sol_pda(end) - T);
t_cyc    = t_sol_pda(mask_cyc);
X_cyc    = X_sol_pda(mask_cyc, :);

% Subplot 1: LV, Aorta, PA pressure
subplot(2, 2, 1);
hold on;
plot(t_cyc, X_cyc(:, pda_idx.P_lv),  '-',  'Color', col_lv,  'LineWidth', 1.5, 'DisplayName', 'P_{LV}');
plot(t_cyc, X_cyc(:, pda_idx.P_ao),  '-',  'Color', col_ao,  'LineWidth', 1.5, 'DisplayName', 'P_{Ao} (proximal)');
plot(t_cyc, X_cyc(:, pda_idx.P_pa),  '--', 'Color', col_pa,  'LineWidth', 1.5, 'DisplayName', 'P_{PA}');
hold off;
xlabel('Time [s]',        'FontSize', font_sz_lb, 'FontName', font_name);
ylabel('Pressure [mmHg]', 'FontSize', font_sz_lb, 'FontName', font_name);
title('Ventricular & Aortic Pressures', 'FontSize', font_sz_ti, 'FontName', font_name);
legend('FontSize', font_sz_ax, 'Location', 'northeast');
grid on; set(gca, 'FontSize', font_sz_ax, 'FontName', font_name, 'Box', 'on');

% Subplot 2: PDA shunt flow
subplot(2, 2, 2);
hold on;
plot(t_cyc, X_cyc(:, pda_idx.Q_shunt_pda), '-', 'Color', col_pda, 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'Zero (no flow)');
hold off;
xlabel('Time [s]',    'FontSize', font_sz_lb, 'FontName', font_name);
ylabel('Flow [mL/s]', 'FontSize', font_sz_lb, 'FontName', font_name);
title(sprintf('PDA Shunt Flow (Mean: %.2f mL/s)', idx_pda.Q_shunt_pda_mean), ...
    'FontSize', font_sz_ti, 'FontName', font_name);
grid on; set(gca, 'FontSize', font_sz_ax, 'FontName', font_name, 'Box', 'on');

% Subplot 3: RA and LA pressures
subplot(2, 2, 3);
hold on;
plot(t_cyc, X_cyc(:, pda_idx.P_ra), '-',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', 'P_{RA}');
plot(t_cyc, X_cyc(:, pda_idx.P_la), '--', 'Color', col_lv,         'LineWidth', 1.5, 'DisplayName', 'P_{LA}');
hold off;
xlabel('Time [s]',        'FontSize', font_sz_lb, 'FontName', font_name);
ylabel('Pressure [mmHg]', 'FontSize', font_sz_lb, 'FontName', font_name);
title('Atrial Pressures', 'FontSize', font_sz_ti, 'FontName', font_name);
legend('FontSize', font_sz_ax, 'Location', 'best');
grid on; set(gca, 'FontSize', font_sz_ax, 'FontName', font_name, 'Box', 'on');

% Subplot 4: Summary text panel
subplot(2, 2, 4); axis off;
summary = {
    sprintf('\\bf PDA Baseline Summary');
    ' ';
    sprintf('HR:        %d bpm',        params_pda.HR_bpm);
    sprintf('R_{PDA}:   %.3f mmHg\cdots/mL', params_pda.R_shunt_pda);
    ' ';
    sprintf('P_{ao} Sys/Dia: %.1f / %.1f mmHg', idx_pda.P_ao_sys, idx_pda.P_ao_dia);
    sprintf('P_{PA} mean:    %.1f mmHg',          idx_pda.P_pa_mean);
    sprintf('SV:             %.2f mL',             idx_pda.SV_lv);
    sprintf('CO:             %.2f L/min',           idx_pda.CO_Lmin);
    sprintf('Qp/Qs:          %.2f',                 idx_pda.Qp_Qs);
    sprintf('EF_{LV}:        %.1f%%',               idx_pda.EF_lv * 100);
    sprintf('SW_{LV}:        %.2f J',               idx_pda.SW_lv_J);
};
text(0.05, 0.95, summary, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', font_sz_ax, 'FontName', font_name, 'Interpreter', 'tex');

%% =========================================================================
%  FIGURE 2: PV LOOP (PDA BASELINE)
% =========================================================================
% Reconstruct E_lv trace for V_lv computation
tn_cyc = mod(t_cyc, T);
En_cyc = zeros(size(tn_cyc));
for k = 1:length(tn_cyc)
    if tn_cyc(k) <= params_pda.Ts1
        En_cyc(k) = 1.55 * (tn_cyc(k) / params_pda.Ts1)^2;
    elseif tn_cyc(k) <= params_pda.Ts2
        alpha = (pi/2)*(tn_cyc(k)-params_pda.Ts1)/(params_pda.Ts2-params_pda.Ts1);
        En_cyc(k) = 1.55 * cos(alpha)^2;
    end
end
E_lv_cyc = (params_pda.Emax_lv - params_pda.Emin_lv)*En_cyc + params_pda.Emin_lv;
V_lv_cyc = X_cyc(:, pda_idx.P_lv) ./ E_lv_cyc + params_pda.V0_lv;

fig2 = figure('Units', 'centimeters', 'Position', [0 0 14 12], 'Color', 'w');
plot(V_lv_cyc, X_cyc(:, pda_idx.P_lv), '-', 'Color', col_lv, 'LineWidth', 2);
xlabel('V_{LV} [mL]',    'FontSize', font_sz_lb, 'FontName', font_name);
ylabel('P_{LV} [mmHg]',  'FontSize', font_sz_lb, 'FontName', font_name);
title(sprintf('LV Pressure-Volume Loop (PDA Baseline) — %s', patient_id), ...
    'FontSize', font_sz_ti, 'FontName', font_name);
grid on; set(gca, 'FontSize', font_sz_ax, 'FontName', font_name, 'Box', 'on');

if ~isempty(coa_scenarios)
    %% =====================================================================
    %  FIGURE 3: TRANS-CoA GRADIENT TRACES (PDA_CoA mode only)
    % =====================================================================
    fig3 = figure('Units', 'centimeters', 'Position', [0 0 22 10], 'Color', 'w', ...
        'Name', 'Trans_CoA_Gradient');
    hold on;
    for s = 1:length(coa_scenarios)
        sc      = coa_scenarios{s};
        sc_idx  = sc.params.idx;
        T_s     = sc.params.T_cardiac;
        mask_s  = sc.t_sol >= (sc.t_sol(end) - T_s);
        t_s     = sc.t_sol(mask_s);
        X_s     = sc.X_sol(mask_s, :);
        dP_s    = X_s(:, sc_idx.P_ao) - X_s(:, sc_idx.P_ao_dist);
        ls      = coa_line_styles{mod(s-1, 3)+1};
        plot(t_s - t_s(1), dP_s, ls, 'LineWidth', 1.5, 'DisplayName', sc.label);
    end
    hold off;
    xlabel('Time in cycle [s]',     'FontSize', font_sz_lb, 'FontName', font_name);
    ylabel('\DeltaP_{CoA} [mmHg]',  'FontSize', font_sz_lb, 'FontName', font_name);
    title(sprintf('Trans-CoA Pressure Gradient — Patient %s', patient_id), ...
        'FontSize', font_sz_ti, 'FontName', font_name);
    legend('FontSize', font_sz_ax, 'Location', 'northeast');
    grid on; set(gca, 'FontSize', font_sz_ax, 'FontName', font_name, 'Box', 'on');

    %% =====================================================================
    %  FIGURE 4: BAR CHART COMPARISON — PDA vs CoA SCENARIOS (PDA_CoA only)
    % =====================================================================
    all_scenarios  = [{struct('indices', results_pda, 'label', 'PDA Only')}, coa_scenarios];
    n_sc           = length(all_scenarios);
    labels         = cell(1, n_sc);
    MAP_vals       = zeros(1, n_sc);
    CO_vals        = zeros(1, n_sc);
    SW_vals        = zeros(1, n_sc);
    dP_coa_vals    = NaN(1, n_sc);   % NaN so PDA-only baseline bar is absent
    Qp_Qs_vals     = zeros(1, n_sc);

    for s = 1:n_sc
        sc_i           = all_scenarios{s};
        labels{s}      = sc_i.label;
        MAP_vals(s)    = sc_i.indices.model.P_ao_mean;
        CO_vals(s)     = sc_i.indices.model.CO_Lmin;
        SW_vals(s)     = sc_i.indices.model.SW_lv_J * 1000;   % → mJ for display
        dP_coa_vals(s) = sc_i.indices.model.DeltaP_coa_mean;  % NaN for PDA-only row
        Qp_Qs_vals(s)  = sc_i.indices.model.Qp_Qs;
    end

    fig4 = figure('Units', 'centimeters', 'Position', [0 0 24 16], 'Color', 'w', ...
        'Name', 'Scenario_Comparison');
    sgtitle(sprintf('Scenario Comparison — Patient %s', patient_id), ...
        'FontSize', font_sz_ti, 'FontName', font_name, 'FontWeight', 'bold');

    metrics = {MAP_vals, CO_vals, SW_vals, dP_coa_vals, Qp_Qs_vals};
    ylabels = {'P_{ao,mean} [mmHg]', 'CO [L/min]', 'SW_{LV} [mJ]', ...
               '\DeltaP_{CoA} [mmHg]', 'Qp/Qs [-]'};
    titles  = {'Mean Arterial Pressure', 'Cardiac Output', ...
               'LV Stroke Work', 'Mean Trans-CoA Gradient', 'Qp/Qs Ratio'};

    for m_idx = 1:5
        subplot(2, 3, m_idx);
        bar(1:n_sc, metrics{m_idx}, 'FaceColor', col_ao, 'EdgeColor', 'k');
        set(gca, 'XTick', 1:n_sc, 'XTickLabel', labels, ...
            'FontSize', font_sz_ax, 'FontName', font_name, 'Box', 'on');
        xtickangle(20);
        ylabel(ylabels{m_idx}, 'FontSize', font_sz_lb, 'FontName', font_name);
        title(titles{m_idx},   'FontSize', font_sz_ti, 'FontName', font_name);
        grid on;
    end

else
    %% =====================================================================
    %  FIGURE 3 (PDA-only): PRIMARY OUTPUTS COMPARISON — MAP and SV
    %  Only generated for PDA-only patients (no CoA scenarios).
    %  Title: "Primary LPM Outputs: MAP and SV"
    % =====================================================================
    clin_vals  = results_pda.clinical;
    m_pda_out  = results_pda.model;

    fig3_pda = figure('Units', 'centimeters', 'Position', [0 0 18 10], 'Color', 'w', ...
        'Name', 'Primary_LPM_Outputs_MAP_SV');
    sgtitle(sprintf('Primary LPM Outputs: MAP and SV — Patient %s', patient_id), ...
        'FontSize', font_sz_ti, 'FontName', font_name, 'FontWeight', 'bold');

    % MAP: clinical vs simulated
    subplot(1, 2, 1);
    bar_map = bar([clin_vals.P_ao_mean_mmHg, m_pda_out.P_ao_mean], ...
        'FaceColor', col_ao, 'EdgeColor', 'k');  %#ok<NASGU>
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Clinical', 'Simulated'}, ...
        'FontSize', font_sz_ax, 'FontName', font_name, 'Box', 'on');
    ylabel('MAP [mmHg]', 'FontSize', font_sz_lb, 'FontName', font_name);
    title('Mean Arterial Pressure', 'FontSize', font_sz_ti, 'FontName', font_name);
    grid on;
    map_pct_err = 100 * abs(m_pda_out.P_ao_mean - clin_vals.P_ao_mean_mmHg) ...
                      / clin_vals.P_ao_mean_mmHg;
    ylim_map = ylim;
    text(1.5, ylim_map(1) + 0.92*(ylim_map(2)-ylim_map(1)), ...
        sprintf('Err: %.1f%%', map_pct_err), ...
        'HorizontalAlignment', 'center', 'FontSize', font_sz_ax, 'FontName', font_name);

    % SV: clinical vs simulated
    subplot(1, 2, 2);
    bar_sv = bar([clin_vals.SV_mL, m_pda_out.SV_lv], ...
        'FaceColor', col_lv, 'EdgeColor', 'k');  %#ok<NASGU>
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Clinical', 'Simulated'}, ...
        'FontSize', font_sz_ax, 'FontName', font_name, 'Box', 'on');
    ylabel('SV [mL]', 'FontSize', font_sz_lb, 'FontName', font_name);
    title('Stroke Volume', 'FontSize', font_sz_ti, 'FontName', font_name);
    grid on;
    sv_pct_err = 100 * abs(m_pda_out.SV_lv - clin_vals.SV_mL) / clin_vals.SV_mL;
    ylim_sv = ylim;
    text(1.5, ylim_sv(1) + 0.92*(ylim_sv(2)-ylim_sv(1)), ...
        sprintf('Err: %.1f%%', sv_pct_err), ...
        'HorizontalAlignment', 'center', 'FontSize', font_sz_ax, 'FontName', font_name);

end

end
