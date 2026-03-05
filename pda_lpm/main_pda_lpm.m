%% MAIN_PDA_LPM
% =================================================================
% PDA LUMPED PARAMETER MODEL — VIRTUAL CoA SEVERITY FRAMEWORK
% Based on: Ortiz-Rangel et al. (2022), Biomed Signal Process Control 71:103151
%
% PURPOSE:
%   1. Calibrate a patient-specific LPM to PDA clinical data
%   2. Simulate PDA-only cardiovascular physiology (baseline)
%   3. Introduce virtual CoA stenosis (50%, 75%, 90%) and simulate
%      resulting haemodynamic changes
%   4. Compare pressure gradients, flow redistribution, and LV workload
%
% USAGE: Run this script from the pda_lpm/ directory.
%
% RULES:
%   - No physics computed in this file (Guardrail §2 and §13.4)
%   - No hardcoded patient values (Guardrail §9.2)
%   - All parameters passed via structs (Guardrail §5.1)
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%   [2] Keshavarz-Motamed et al. (2011). J Biomech 44:2817–2825.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% =================================================================

clear; clc; close all;

%% Add all subdirectory paths
addpath('config', 'models', 'solvers', 'utils', 'tests');

fprintf('=================================================================\n');
fprintf('   PDA LPM — VIRTUAL CoA SEVERITY FRAMEWORK\n');
fprintf('   Ortiz-Rangel et al. (2022), Biomed Signal Process Control\n');
fprintf('=================================================================\n\n');

%% =========================================================================
%  STEP 1 — LOAD PATIENT DATA
% =========================================================================
fprintf('STEP 1: Loading patient data...\n');

csv_path = fullfile('config', 'patient_data.csv');
[clinical, ~] = load_patient_data(csv_path);

%% =========================================================================
%  STEP 2 — BUILD PATIENT-SPECIFIC PARAMETERS
% =========================================================================
fprintf('STEP 2: Building patient-specific LPM parameters...\n');

params_default  = default_parameters();
params_pda      = build_patient_params(clinical, params_default);

%% =========================================================================
%  STEP 3 — SIMULATE PDA-ONLY BASELINE
% =========================================================================
fprintf('STEP 3: Simulating PDA-only cardiovascular model (baseline)...\n');

rhs_pda         = @(t, X) system_rhs_pda(t, X, params_pda);
n_warmup        = 8;    % Warm-up cycles before reporting (Guardrail §8.3)
n_report        = 2;    % Steady-state cycles to analyse

[t_sol_pda, X_sol_pda] = integrate_system(rhs_pda, params_pda, n_warmup, n_report);

results_pda = compute_clinical_indices(t_sol_pda, X_sol_pda, params_pda, ...
                                       clinical, 'PDA Only (Baseline)');

%% =========================================================================
%  STEP 4 — SELECT CoA STENOSIS SCENARIOS
% =========================================================================
fprintf('STEP 4: Select virtual CoA stenosis severities to simulate.\n');
fprintf('  Available: 50%% (mild-moderate), 75%% (moderate-severe), 90%% (critical)\n');
fprintf('  Enter severities as a vector, e.g., [50 75 90]:\n');

stenosis_input = input('  CoA stenosis percentages: ');
if isempty(stenosis_input)
    stenosis_input = [50, 75, 90];
    fprintf('  (Using default: [50 75 90])\n');
end

%% =========================================================================
%  STEP 5 — SIMULATE EACH CoA SCENARIO
% =========================================================================
fprintf('\nSTEP 5: Simulating virtual CoA scenarios...\n\n');

coa_scenarios = cell(1, length(stenosis_input));

for s_idx = 1:length(stenosis_input)
    s_pct   = stenosis_input(s_idx);
    s_label = sprintf('CoA %.0f%% Stenosis', s_pct);

    fprintf('--- Scenario: %s ---\n', s_label);

    % Build CoA parameter set (extends PDA params)
    params_coa_s = build_coa_params(params_pda, clinical, s_pct);

    % Simulate with PDA + CoA RHS
    rhs_coa_s    = @(t, X) system_rhs_pda_coa(t, X, params_coa_s);
    [t_s, X_s]   = integrate_system(rhs_coa_s, params_coa_s, n_warmup, n_report);

    % Compute clinical indices for this scenario
    idx_s = compute_clinical_indices(t_s, X_s, params_coa_s, clinical, s_label);

    % Store scenario results
    coa_scenarios{s_idx} = struct( ...
        't_sol',   t_s,        ...
        'X_sol',   X_s,        ...
        'indices', idx_s,      ...
        'params',  params_coa_s, ...
        'label',   s_label);
end

%% =========================================================================
%  STEP 6 — GENERATE PUBLICATION-READY FIGURES
% =========================================================================
fprintf('STEP 6: Generating publication-ready figures...\n');

plot_results(t_sol_pda, X_sol_pda, results_pda, coa_scenarios, clinical.patient_id);

%% =========================================================================
%  STEP 7 — PRINT CONSOLIDATED COMPARISON TABLE
% =========================================================================
fprintf('\nSTEP 7: Consolidated scenario comparison\n');
fprintf('=================================================================\n');
fprintf('%-22s %8s %8s %8s %10s %8s %8s\n', ...
    'Scenario', 'MAP', 'CO', 'SV', 'DeltaP_CoA', 'EF_LV', 'Qp/Qs');
fprintf('%-22s %8s %8s %8s %10s %8s %8s\n', ...
    '', 'mmHg', 'L/min', 'mL', 'mmHg', '%', '-');
fprintf('%s\n', repmat('-', 1, 76));

% PDA-only row
m = results_pda.model;
fprintf('%-22s %8.1f %8.2f %8.2f %10.1f %8.1f %8.2f\n', ...
    'PDA Only', m.P_ao_mean, m.CO_Lmin, m.SV_lv, ...
    m.DeltaP_coa_mean, m.EF_lv*100, m.Qp_Qs);

% CoA scenario rows
for s_idx = 1:length(coa_scenarios)
    sc  = coa_scenarios{s_idx};
    mc  = sc.indices.model;
    fprintf('%-22s %8.1f %8.2f %8.2f %10.1f %8.1f %8.2f\n', ...
        sc.label, mc.P_ao_mean, mc.CO_Lmin, mc.SV_lv, ...
        mc.DeltaP_coa_mean, mc.EF_lv*100, mc.Qp_Qs);
end

fprintf('%s\n', repmat('=', 1, 76));
fprintf('  Clinical reference (measured):  MAP = %.1f mmHg | SV = %.2f mL\n', ...
    clinical.P_ao_mean_mmHg, clinical.SV_mL);
fprintf('=================================================================\n\n');
fprintf('Simulation complete.\n');
fprintf('Figures: see on-screen.\n');
fprintf('To export figures as PDF: set MATLAB export to vector PDF.\n');
