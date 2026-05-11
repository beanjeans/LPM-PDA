%% MAIN_PDA_LPM
% =================================================================
% PDA LUMPED PARAMETER MODEL — VIRTUAL CoA SEVERITY FRAMEWORK
% Based on: Ortiz-Rangel et al. (2022), Biomed Signal Process Control 71:103151
%
% PURPOSE:
%   1. Calibrate a patient-specific LPM to PDA clinical data
%   2. Simulate PDA-only cardiovascular physiology (baseline)
%   3. Introduce virtual CoA (stenosis × length combinations) and simulate
%      resulting haemodynamic changes
%   4. Report CoA severity based on SIMULATED PRESSURE GRADIENT (not
%      stenosis_pct directly), following ESC guideline thresholds
%
% KEY DESIGN PRINCIPLE:
%   stenosis_pct  → drives CoA geometry and pressure gradient
%   ΔP_CoA        → drives predicted_CoA_severity classification
%   PDA           → retained as haemodynamic modifier (may mask ΔP_CoA)
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
%   [3] Baumgartner et al. (2010). Eur Heart J 31(19):2369–2417.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  3.0  — allometric scaling (Pennati & Fumero 2000); variable CoA length
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

% Convert clinical weight from grams to kg for allometric scaling
BW_neo_kg = clinical.weight_g / 1000;   % [g] → [kg]
% Gestational age: read from CSV if available, otherwise use default
if isfield(clinical, 'GA_weeks')
    GA_neo_weeks = clinical.GA_weeks;
else
    GA_neo_weeks = 29;   % default preterm reference [weeks]
end

% Build allometrically-scaled default parameters for this patient's body weight
params_default  = default_parameters(BW_neo_kg, GA_neo_weeks);
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
%  STEP 4 — SELECT CoA STENOSIS SCENARIOS AND LENGTHS
% =========================================================================
fprintf('STEP 4: Select virtual CoA stenosis severities and lengths to simulate.\n');
fprintf('  Stenosis:  e.g. [50 75 90] (%% area stenosis)\n');
fprintf('  Length:    e.g. [3 8]  (mm) — discrete(<5mm) or long-segment(>=5mm)\n\n');

stenosis_input = input('  CoA stenosis percentages [default: 50 75 90]: ');
if isempty(stenosis_input)
    stenosis_input = [50, 75, 90];
    fprintf('  (Using default stenosis: [50 75 90]%%)\n');
end

length_input = input('  CoA length(s) in mm    [default: 3 8]:  ');
if isempty(length_input)
    length_input = [3, 8];   % 3 mm = discrete; 8 mm = long-segment
    fprintf('  (Using default lengths: [3 8] mm)\n');
end

%% =========================================================================
%  STEP 5 — SIMULATE EACH CoA SCENARIO (stenosis × length combinations)
% =========================================================================
fprintf('\nSTEP 5: Simulating virtual CoA scenarios...\n\n');

n_scenarios  = length(stenosis_input) * length(length_input);
coa_scenarios = cell(1, n_scenarios);
sc_count = 0;

for len_idx = 1:length(length_input)
    l_mm  = length_input(len_idx);
    for s_idx = 1:length(stenosis_input)
        s_pct   = stenosis_input(s_idx);
        sc_count = sc_count + 1;

        if l_mm < 5
            len_label = sprintf('%.0fmm-discrete', l_mm);
        else
            len_label = sprintf('%.0fmm-longseg', l_mm);
        end
        s_label = sprintf('CoA %.0f%% | %s', s_pct, len_label);

        fprintf('--- Scenario %d/%d: %s ---\n', sc_count, n_scenarios, s_label);

        % Build CoA parameter set (extends PDA params)
        % Pass coa_length_mm as 4th argument (new in v2.0)
        params_coa_s = build_coa_params(params_pda, clinical, s_pct, l_mm);

        % Simulate with PDA + CoA RHS
        rhs_coa_s    = @(t, X) system_rhs_pda_coa(t, X, params_coa_s);
        [t_s, X_s]   = integrate_system(rhs_coa_s, params_coa_s, n_warmup, n_report);

        % Compute clinical indices (severity classification inside)
        idx_s = compute_clinical_indices(t_s, X_s, params_coa_s, clinical, s_label);

        % Store scenario results
        coa_scenarios{sc_count} = struct( ...
            't_sol',   t_s,        ...
            'X_sol',   X_s,        ...
            'indices', idx_s,      ...
            'params',  params_coa_s, ...
            'label',   s_label);
    end
end

%% =========================================================================
%  STEP 6 — GENERATE PUBLICATION-READY FIGURES
% =========================================================================
fprintf('STEP 6: Generating publication-ready figures...\n');

plot_results(t_sol_pda, X_sol_pda, results_pda, coa_scenarios, clinical.patient_id);

%% =========================================================================
%  STEP 7 — PRINT CONSOLIDATED CoA SEVERITY TABLE
%  PRIMARY OUTPUT: predicted_CoA_severity (based on simulated ΔP_CoA)
%  stenosis_pct drives geometry only — NOT reported as severity directly.
% =========================================================================
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
fprintf('=========================================================================================\n\n');

%% =========================================================================
%  STEP 8 — PRINT HAEMODYNAMIC COMPARISON TABLE (MAP / CO / SV)
% =========================================================================
fprintf('STEP 8: Haemodynamic comparison table\n');
fprintf('%s\n', repmat('-', 1, 76));
fprintf('%-26s %8s %8s %8s %10s %8s %8s\n', ...
    'Scenario', 'MAP', 'CO', 'SV', 'DeltaP_CoA', 'EF_LV', 'Qp/Qs');
fprintf('%-26s %8s %8s %8s %10s %8s %8s\n', ...
    '', 'mmHg', 'L/min', 'mL', 'mmHg', '%', '-');
fprintf('%s\n', repmat('-', 1, 76));

% PDA-only row
m = results_pda.model;
fprintf('%-26s %8.1f %8.2f %8.2f %10.1f %8.1f %8.2f\n', ...
    'PDA Only', m.P_ao_mean, m.CO_Lmin, m.SV_lv, ...
    m.DeltaP_coa_mean, m.EF_lv*100, m.Qp_Qs);

% CoA scenario rows
for s_idx = 1:length(coa_scenarios)
    sc  = coa_scenarios{s_idx};
    mc  = sc.indices.model;
    fprintf('%-26s %8.1f %8.2f %8.2f %10.1f %8.1f %8.2f\n', ...
        sc.label, mc.P_ao_mean, mc.CO_Lmin, mc.SV_lv, ...
        mc.DeltaP_coa_mean, mc.EF_lv*100, mc.Qp_Qs);
end

fprintf('%s\n', repmat('=', 1, 76));
fprintf('  Clinical reference:  MAP = %.1f mmHg | SV = %.2f mL\n', ...
    clinical.P_ao_mean_mmHg, clinical.SV_mL);
fprintf('%s\n', repmat('=', 1, 76));
fprintf('\nSimulation complete.\n');
fprintf('Figures: see on-screen.\n');
fprintf('To export figures as PDF: set MATLAB export to vector PDF.\n');

