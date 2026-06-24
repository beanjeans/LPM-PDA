%% MAIN_PDA_LPM
% =================================================================
% PDA LUMPED PARAMETER MODEL — FINAL SIMULATION
% Based on: Ortiz-Rangel et al. (2022), Biomed Signal Process Control 71:103151
%
% PURPOSE:
%   1. Load patient data and auto-detect disease mode (PDA_only vs PDA_CoA).
%   2. Build patient-specific LPM parameters; optionally apply optimised
%      parameters from results/optimization/optimization_workspace.mat.
%   3. Simulate PDA-only cardiovascular physiology (baseline) using
%      system_rhs_pda.m.
%   4. For PDA-CoA patients only: simulate virtual CoA scenarios
%      (stenosis × length combinations) using system_rhs_pda_coa.m, then
%      classify CoA severity from the simulated pressure gradient.
%
% DISEASE MODE DETECTION (automatic — no user input required):
%   PDA_only — dP_coa_mmHg <= 0  AND  D_coa_mm <= 0  AND  v_coa_ms <= 0
%              RHS: system_rhs_pda.m   CoA module: NOT activated.
%              Final outputs: MAP and SV only.
%   PDA_CoA  — any CoA measurement present in the clinical record.
%              RHS: system_rhs_pda_coa.m for each CoA scenario.
%              Final outputs: full haemodynamic + CoA severity table.
%
% OUTPUT SELECTION (P01 / P02 — PDA_only):
%   - Final table shows MAP and SV only.
%   - CoA gradient, severity, and stenosis_pct are NOT reported.
%   - stenosis_pct from optimization_workspace.mat is silently ignored.
%   - CoA plots are not generated.
%
% OUTPUT SELECTION (P03 — PDA_CoA):
%   - Full CoA scenario severity table and haemodynamic comparison.
%   - Existing PDA-CoA workflow preserved unchanged.
%
% USAGE: Run this script from the pda_lpm/ directory.
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%   [2] Keshavarz-Motamed et al. (2011). J Biomech 44:2817–2825.
%   [3] Baumgartner et al. (2010). Eur Heart J 31(19):2369–2417.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  4.0  — disease-mode routing; PDA-only / PDA-CoA conditional pipeline
% =================================================================

clear; clc; close all;

%% Add all subdirectory paths
addpath('config', 'models', 'solvers', 'utils', 'tests', 'optimization');

fprintf('=================================================================\n');
fprintf('   PDA LPM — FINAL SIMULATION\n');
fprintf('   Ortiz-Rangel et al. (2022), Biomed Signal Process Control\n');
fprintf('=================================================================\n\n');

%% =========================================================================
%  STEP 1 — LOAD PATIENT DATA
% =========================================================================
fprintf('STEP 1: Loading patient data...\n');

csv_path = fullfile('config', 'patient_data.csv');
[clinical, ~] = load_patient_data(csv_path);

%% =========================================================================
%  STEP 2 — DETECT DISEASE MODE
%  Criteria for PDA_only: no CoA pressure gradient, no CoA diameter,
%  and no CoA Doppler velocity recorded (all <= 0 or zero in CSV).
%  Patients with any positive CoA measurement are routed to PDA_CoA mode.
% =========================================================================
fprintf('STEP 2: Detecting disease mode...\n');

is_pda_only = (clinical.dP_coa_mmHg <= 0) && ...
              (clinical.D_coa_mm     <= 0) && ...
              (clinical.v_coa_ms     <= 0);

if is_pda_only
    disease_mode = 'PDA_only';
    fprintf('  Disease mode: PDA_only\n');
    fprintf('  PDA-only patient detected. Applying optimized PDA/systemic parameters only: Emax_lv and R_systemic.\n');
    fprintf('  CoA geometry parameters are ignored for this patient.\n\n');
else
    disease_mode = 'PDA_CoA';
    fprintf('  Disease mode: PDA_CoA (CoA measurements detected — virtual CoA simulation enabled).\n\n');
end

%% =========================================================================
%  STEP 3 — BUILD PATIENT-SPECIFIC PARAMETERS
% =========================================================================
fprintf('STEP 3: Building patient-specific LPM parameters...\n');

BW_neo_kg      = clinical.weight_g / 1000;   % [g] → [kg]
params_default = default_parameters(BW_neo_kg);
params_pda     = build_patient_params(clinical, params_default);

%% =========================================================================
%  STEP 3b — APPLY OPTIMISED PARAMETERS (if available)
%  For PDA_only patients: apply only Emax_lv and R_systemic.
%  stenosis_pct and coa_length_mm are silently ignored — they affect the
%  CoA module only, which is not activated for PDA-only patients.
% =========================================================================
opt_workspace = fullfile('results', 'optimization', 'optimization_workspace.mat');

if exist(opt_workspace, 'file')
    fprintf('  Loading optimised parameters from: %s\n', opt_workspace);
    ws_opt = load(opt_workspace);

    if isfield(ws_opt, 'x_opt') && isfield(ws_opt, 'opt_config')
        x_opt_loaded         = ws_opt.x_opt;
        opt_cfg              = ws_opt.opt_config;
        opt_cfg.fixed_params = params_pda;   % use freshly built patient params as base

        [params_opt_raw, ~, ~] = apply_optimized_params(x_opt_loaded, opt_cfg);

        if strcmp(disease_mode, 'PDA_only')
            % PDA-only: copy only the physically relevant optimised params.
            % stenosis_pct and coa_length_mm are NOT applied — CoA module is off.
            params_pda.Emax_lv    = params_opt_raw.Emax_lv;
            params_pda.Emin_lv    = params_opt_raw.Emin_lv;
            params_pda.Emax_rv    = params_opt_raw.Emax_rv;
            params_pda.Emin_rv    = params_opt_raw.Emin_rv;
            params_pda.R_systemic = params_opt_raw.R_systemic;
            fprintf('  Applied optimised: Emax_lv = %.4f mmHg/mL,  R_systemic = %.4f mmHg*s/mL\n', ...
                params_pda.Emax_lv, params_pda.R_systemic);
            fprintf('  (stenosis_pct and coa_length_mm ignored — CoA module not active for PDA-only patient)\n\n');
        else
            params_pda = params_opt_raw;
            fprintf('  Applied all optimised parameters for PDA_CoA mode.\n\n');
        end
    else
        fprintf('  WARNING: optimization_workspace.mat missing x_opt or opt_config — using default calibrated params.\n\n');
    end
else
    fprintf('  No optimisation workspace found — using default calibrated parameters.\n\n');
end

%% =========================================================================
%  STEP 4 — SIMULATE PDA-ONLY BASELINE
%  Always uses system_rhs_pda.m (no CoA module) regardless of disease mode.
% =========================================================================
fprintf('STEP 4: Simulating PDA-only cardiovascular model (baseline)...\n');

rhs_pda     = @(t, X) system_rhs_pda(t, X, params_pda);
n_warmup    = 8;    % Warm-up cycles before reporting (transient settling)
n_report    = 2;    % Steady-state cycles to analyse

[t_sol_pda, X_sol_pda] = integrate_system(rhs_pda, params_pda, n_warmup, n_report);

results_pda = compute_clinical_indices(t_sol_pda, X_sol_pda, params_pda, ...
                                       clinical, 'PDA Only (Baseline)');

%% =========================================================================
%  STEP 5 — VIRTUAL CoA SCENARIOS (PDA_CoA patients only)
%  For PDA_only patients this step is entirely skipped.
%  coa_scenarios is always initialised as an empty cell array so that
%  downstream code (plot_results, save loop) is safe regardless of mode.
% =========================================================================
coa_scenarios = {};   % Always initialise as empty cell array

if strcmp(disease_mode, 'PDA_CoA')

    fprintf('STEP 5: Select virtual CoA stenosis severities and lengths to simulate.\n');
    fprintf('  Stenosis:  e.g. [50 75 90] (%% area stenosis)\n');
    fprintf('  Length:    e.g. [3 8]  (mm) — discrete(<5mm) or long-segment(>=5mm)\n\n');

    stenosis_input = input('  CoA stenosis percentages [default: 50 75 90]: ');
    if isempty(stenosis_input)
        stenosis_input = [50, 75, 90];
        fprintf('  (Using default stenosis: [50 75 90]%%)\n');
    end

    length_input = input('  CoA length(s) in mm    [default: 3 8]:  ');
    if isempty(length_input)
        length_input = [3, 8];
        fprintf('  (Using default lengths: [3 8] mm)\n');
    end

    fprintf('\nSimulating virtual CoA scenarios...\n\n');

    n_scenarios = length(stenosis_input) * length(length_input);
    sc_count    = 0;

    for len_idx = 1:length(length_input)
        l_mm = length_input(len_idx);
        for s_idx = 1:length(stenosis_input)
            s_pct    = stenosis_input(s_idx);
            sc_count = sc_count + 1;

            if l_mm < 5
                len_label = sprintf('%.0fmm-discrete', l_mm);
            else
                len_label = sprintf('%.0fmm-longseg', l_mm);
            end
            s_label = sprintf('CoA %.0f%% | %s', s_pct, len_label);

            fprintf('--- Scenario %d/%d: %s ---\n', sc_count, n_scenarios, s_label);

            % Build CoA parameter set and simulate with PDA+CoA RHS
            params_coa_s = build_coa_params(params_pda, clinical, s_pct, l_mm);
            rhs_coa_s    = @(t, X) system_rhs_pda_coa(t, X, params_coa_s);
            [t_s, X_s]   = integrate_system(rhs_coa_s, params_coa_s, n_warmup, n_report);
            idx_s        = compute_clinical_indices(t_s, X_s, params_coa_s, clinical, s_label);

            coa_scenarios{sc_count} = struct( ...
                't_sol',   t_s,          ...
                'X_sol',   X_s,          ...
                'indices', idx_s,        ...
                'params',  params_coa_s, ...
                'label',   s_label);
        end
    end

else
    fprintf('STEP 5: PDA-only patient — skipping virtual CoA scenario simulation.\n');
    fprintf('  No CoA module is applied for this patient.\n\n');
end

%% =========================================================================
%  STEP 6 — GENERATE FIGURES
%  plot_results receives params_pda as 6th argument so it can safely
%  access state indices even when coa_scenarios is empty.
% =========================================================================
fprintf('STEP 6: Generating figures...\n');

plot_results(t_sol_pda, X_sol_pda, results_pda, coa_scenarios, clinical.patient_id, params_pda);

%% =========================================================================
%  STEP 7 — CoA SEVERITY TABLE (PDA_CoA patients only)
%  For PDA_only patients this step is skipped entirely.
%  CoA gradient is never shown as 0 mmHg — it is simply not reported.
% =========================================================================
if strcmp(disease_mode, 'PDA_CoA') && ~isempty(coa_scenarios)

    fprintf('\nSTEP 7: CoA Clinical Severity Summary\n');
    fprintf('=========================================================================================\n');
    fprintf('%-26s %6s %7s %10s %11s %10s %8s %s\n', ...
        'Scenario', 'S_pct', 'L_mm', 'Category', 'dP_peak', 'dP_mean_s', 'Qcoa/Qt', 'Severity');
    fprintf('%-26s %6s %7s %10s %11s %10s %8s %s\n', ...
        '', '%', 'mm', '', 'mmHg', 'mmHg', '-', '');
    fprintf('%s\n', repmat('-', 1, 89));

    for s_idx = 1:length(coa_scenarios)
        sc  = coa_scenarios{s_idx};
        mc  = sc.indices.model;
        fprintf('%-26s %6.0f %7.1f %10s %11.1f %10.1f %8.2f %s\n', ...
            sc.label, ...
            mc.stenosis_pct, ...
            mc.coa_length_mm, ...
            mc.coa_length_category, ...
            mc.DeltaP_coa_peak, ...
            mc.DeltaP_coa_mean_sys, ...
            mc.Q_coa_fraction, ...
            upper(mc.predicted_CoA_severity));
    end

    fprintf('%s\n', repmat('=', 1, 89));
    fprintf('  Severity thresholds (ESC Ref [3]):  Mild <20 mmHg | Moderate 20-40 mmHg | Severe >40 mmHg\n');
    fprintf('  NOTE: stenosis_%% drives geometry only. Severity is reported via simulated dP_CoA.\n');
    fprintf('  NOTE: PDA (if present) may reduce dP_CoA and mask true CoA severity.\n');
    fprintf('%s\n\n', repmat('=', 1, 89));

else
    fprintf('STEP 7: PDA-only patient — no CoA severity classification applied.\n');
    fprintf('  No CoA module is applied for this patient.\n\n');
end

%% =========================================================================
%  STEP 8 — PRIMARY OUTPUT TABLE
%  PDA_only patients: MAP and SV only (direct optimisation targets).
%  PDA_CoA patients:  full haemodynamic comparison including CoA columns.
% =========================================================================
m = results_pda.model;

if strcmp(disease_mode, 'PDA_only')

    fprintf('STEP 8: PRIMARY LPM OUTPUTS — PDA-ONLY PATIENT\n');
    fprintf('%s\n', repmat('-', 1, 72));
    fprintf('%-12s  %14s  %14s  %10s  %14s\n', ...
        'Target', 'Clinical', 'Simulated', 'Error', 'Percent Error');
    fprintf('%-12s  %14s  %14s  %10s  %14s\n', ...
        '', '(reference)', '(model)', '', '(%)');
    fprintf('%s\n', repmat('-', 1, 72));

    map_err = m.P_ao_mean - clinical.P_ao_mean_mmHg;
    sv_err  = m.SV_lv     - clinical.SV_mL;

    fprintf('%-12s  %14.3f  %14.3f  %10.3f  %13.2f%%\n', ...
        'MAP (mmHg)', clinical.P_ao_mean_mmHg, m.P_ao_mean, map_err, ...
        100 * abs(map_err) / clinical.P_ao_mean_mmHg);
    fprintf('%-12s  %14.3f  %14.3f  %10.3f  %13.2f%%\n', ...
        'SV (mL)', clinical.SV_mL, m.SV_lv, sv_err, ...
        100 * abs(sv_err) / clinical.SV_mL);

    fprintf('%s\n', repmat('=', 1, 72));
    fprintf('  PDA-only final LPM simulation. Primary outputs: MAP and SV.\n');
    fprintf('  No CoA module is applied for this patient.\n\n');

else

    fprintf('STEP 8: Haemodynamic comparison table\n');
    fprintf('%s\n', repmat('-', 1, 76));
    fprintf('%-26s %8s %8s %8s %10s %8s %8s\n', ...
        'Scenario', 'MAP', 'CO', 'SV', 'DeltaP_CoA', 'EF_LV', 'Qp/Qs');
    fprintf('%-26s %8s %8s %8s %10s %8s %8s\n', ...
        '', 'mmHg', 'L/min', 'mL', 'mmHg', '%', '-');
    fprintf('%s\n', repmat('-', 1, 76));

    % PDA-only baseline row — DeltaP_CoA is N/A (no CoA module in baseline)
    fprintf('%-26s %8.1f %8.2f %8.2f %10s %8.1f %8.2f\n', ...
        'PDA Only', m.P_ao_mean, m.CO_Lmin, m.SV_lv, ...
        'N/A', m.EF_lv * 100, m.Qp_Qs);

    for s_idx = 1:length(coa_scenarios)
        sc  = coa_scenarios{s_idx};
        mc  = sc.indices.model;
        fprintf('%-26s %8.1f %8.2f %8.2f %10.1f %8.1f %8.2f\n', ...
            sc.label, mc.P_ao_mean, mc.CO_Lmin, mc.SV_lv, ...
            mc.DeltaP_coa_mean, mc.EF_lv * 100, mc.Qp_Qs);
    end

    fprintf('%s\n', repmat('=', 1, 76));
    fprintf('  Clinical reference:  MAP = %.1f mmHg | SV = %.2f mL\n', ...
        clinical.P_ao_mean_mmHg, clinical.SV_mL);
    fprintf('%s\n', repmat('=', 1, 76));

end

%% =========================================================================
%  STEP 9 — SAVE RESULTS TO FILE
% =========================================================================
fprintf('STEP 9: Saving LPM results...\n');

results_dir_lpm = fullfile('results', 'lpm');
if ~exist(results_dir_lpm, 'dir')
    mkdir(results_dir_lpm);
end

patient_id_safe = strrep(clinical.patient_id, ' ', '_');

if strcmp(disease_mode, 'PDA_only')

    % --- 9a. Save PDA-only workspace ---
    mat_path_lpm = fullfile(results_dir_lpm, sprintf('lpm_results_%s.mat', patient_id_safe));
    save(mat_path_lpm, 't_sol_pda', 'X_sol_pda', 'results_pda', 'params_pda', 'clinical');
    fprintf('  Saved: %s\n', mat_path_lpm);

    % --- 9b. Save primary outputs CSV (MAP and SV only) ---
    m_out      = results_pda.model;
    map_err_o  = m_out.P_ao_mean - clinical.P_ao_mean_mmHg;
    sv_err_o   = m_out.SV_lv     - clinical.SV_mL;

    rows_out = {
        clinical.patient_id, 'MAP', clinical.P_ao_mean_mmHg, m_out.P_ao_mean, ...
            map_err_o, 100 * abs(map_err_o) / clinical.P_ao_mean_mmHg, 'mmHg';
        clinical.patient_id, 'SV',  clinical.SV_mL,          m_out.SV_lv, ...
            sv_err_o,  100 * abs(sv_err_o)  / clinical.SV_mL,          'mL'
    };

    T_lpm = cell2table(rows_out, 'VariableNames', { ...
        'PatientID', 'Metric', 'ClinicalTarget', 'Simulated', ...
        'SignedError', 'PercentError_pct', 'Unit'});

    csv_path_lpm = fullfile(results_dir_lpm, ...
        sprintf('lpm_pda_primary_outputs_%s.csv', patient_id_safe));
    writetable(T_lpm, csv_path_lpm);
    fprintf('  Saved: %s\n', csv_path_lpm);

else

    % --- 9a. Save full workspace (waveforms + indices for all scenarios) ---
    mat_path_lpm = fullfile(results_dir_lpm, sprintf('lpm_results_%s.mat', patient_id_safe));
    save(mat_path_lpm, ...
        't_sol_pda', 'X_sol_pda', 'results_pda', ...
        'coa_scenarios', 'params_pda', 'clinical');
    fprintf('  Saved: %s\n', mat_path_lpm);

    % --- 9b. Save CoA severity summary as CSV ---
    rows_cell = {};
    for s_idx = 1:length(coa_scenarios)
        sc = coa_scenarios{s_idx};
        mc = sc.indices.model;
        rows_cell{end+1, 1} = clinical.patient_id;  %#ok<SAGROW>
        rows_cell{end, 2}   = sc.label;
        rows_cell{end, 3}   = mc.stenosis_pct;
        rows_cell{end, 4}   = mc.coa_length_mm;
        rows_cell{end, 5}   = mc.coa_length_category;
        rows_cell{end, 6}   = mc.DeltaP_coa_peak;
        rows_cell{end, 7}   = mc.DeltaP_coa_mean_sys;
        rows_cell{end, 8}   = mc.Q_coa_fraction;
        rows_cell{end, 9}   = mc.predicted_CoA_severity;
        rows_cell{end, 10}  = mc.P_ao_mean;
        rows_cell{end, 11}  = mc.CO_Lmin;
        rows_cell{end, 12}  = mc.SV_lv;
        rows_cell{end, 13}  = mc.EF_lv * 100;
        rows_cell{end, 14}  = mc.Qp_Qs;
    end

    T_lpm = cell2table(rows_cell, 'VariableNames', { ...
        'PatientID', 'Scenario', 'Stenosis_pct', 'Length_mm', 'Length_category', ...
        'DeltaP_peak_mmHg', 'DeltaP_mean_sys_mmHg', 'Q_coa_fraction', ...
        'Predicted_severity', 'MAP_mmHg', 'CO_Lmin', 'SV_mL', 'EF_pct', 'Qp_Qs'});

    csv_path_lpm = fullfile(results_dir_lpm, sprintf('lpm_coa_summary_%s.csv', patient_id_safe));
    writetable(T_lpm, csv_path_lpm);
    fprintf('  Saved: %s\n', csv_path_lpm);

end

% --- 9c. Export all open figures as PNG ---
fig_handles = findall(0, 'Type', 'figure');
for f = 1:length(fig_handles)
    fig_name = get(fig_handles(f), 'Name');
    if isempty(fig_name)
        fig_name = sprintf('figure_%d', fig_handles(f).Number);
    end
    fig_name_safe = regexprep(fig_name, '[^a-zA-Z0-9_\-]', '_');
    fig_path = fullfile(results_dir_lpm, sprintf('%s_%s.png', patient_id_safe, fig_name_safe));
    exportgraphics(fig_handles(f), fig_path, 'Resolution', 300);
    fprintf('  Saved: %s\n', fig_path);
end

fprintf('\nSimulation complete.\n');
fprintf('Results saved to: %s/\n', results_dir_lpm);
