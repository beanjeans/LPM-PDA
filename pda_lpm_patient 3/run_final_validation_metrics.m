%% RUN_FINAL_VALIDATION_METRICS
% =========================================================================
% DIRECT OPTIMIZATION OUTPUT METRICS — PDA-CoA LPM  (Patient 3)
%
% PURPOSE:
%   Compute and report agreement between the final post-optimisation
%   simulation outputs and the patient's clinical reference values,
%   restricted to the direct clinical targets used in the objective.
%
%   *** SCOPE ***
%   Final validation metrics are computed only for direct optimization
%   outputs: MAP, SV, and dP_CoA_peak.
%   This optimization uses only MAP, SV, and measured dP_CoA_peak as
%   direct clinical targets.  Other model outputs (SBP, DBP, dP_PDA,
%   Q_CoA/Q_total, CO, Qp/Qs, severity) are not included in the
%   optimization or final validation metrics for this configuration.
%
%   *** CALIBRATION NOTE ***
%   These are POST-CALIBRATION agreement metrics, NOT independent external
%   validation.  The same clinical targets used in the objective (MAP, SV,
%   dP_CoA_peak) are evaluated here.  These numbers quantify HOW WELL the
%   model reproduces the training targets, not predictive accuracy on
%   unseen data.
%
% WORKFLOW POSITION:
%   Run AFTER run_lbfgsb_optimization_pda_coa.m (Patient 3 config):
%     1. run_sobol_gsa_pda_coa.m
%     2. run_lbfgsb_optimization_pda_coa.m   ← produces optimization_workspace.mat
%     3. (optional) main_pda_lpm.m
%
%   This script is STANDALONE — it does NOT re-run GSA, optimisation, or
%   simulation.  It only reads the saved workspace and patient CSV.
%
% INPUTS (read from existing files):
%   results/optimization/optimization_workspace.mat
%   config/patient_data.csv
%
% OUTPUTS:
%   results/tables/final_validation_metrics_patient_<ID>.csv
%   results/tables/final_validation_summary_patient_<ID>.csv
%   Console table printed to MATLAB command window.
%
% CSV COLUMNS (per-metric table):
%   Metric, ClinicalTarget, PreOptValue, PostOptValue,
%   SignedError, AbsError, PercentError, Unit, IncludedInObjective
%
% SUMMARY STATISTICS (over MAP, SV, dP_CoA_peak only):
%   MAE, RMSE, MeanPercentError, NMetrics
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%   [2] Baumgartner et al. (2010). Eur Heart J 31(19):2369-2417.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  2.0  — restricted to direct optimization outputs (MAP, SV, dP_CoA_peak)
% =========================================================================

clear; clc;

fprintf('=================================================================\n');
fprintf('   DIRECT OPTIMIZATION OUTPUT METRICS\n');
fprintf('   PDA-CoA Lumped Parameter Cardiovascular Model — Patient 3\n');
fprintf('=================================================================\n\n');
fprintf('  This optimization uses only MAP, SV, and measured dP_CoA_peak\n');
fprintf('  as direct clinical targets.\n');
fprintf('  Other model outputs are not included in the optimization or\n');
fprintf('  final validation metrics for this configuration.\n\n');

% =========================================================================
%  STEP 1 — LOCATE AND LOAD THE OPTIMISATION WORKSPACE
% =========================================================================
fprintf('STEP 1: Loading optimisation workspace...\n');

mat_path = fullfile('results', 'optimization', 'optimization_workspace.mat');

if ~exist(mat_path, 'file')
    fprintf('\n  ERROR: Workspace file not found:\n');
    fprintf('    %s\n\n', mat_path);
    fprintf('  Please run run_lbfgsb_optimization_pda_coa.m first,\n');
    fprintf('  which saves optimization_workspace.mat to results/optimization/.\n\n');
    return;
end

ws = load(mat_path);
fprintf('  Loaded: %s\n\n', mat_path);

% Verify required variables exist
required_ws_vars = {'opt_outputs', 'clinical', 'x_opt', 'J_final'};
missing_ws = {};
for k = 1:length(required_ws_vars)
    if ~isfield(ws, required_ws_vars{k})
        missing_ws{end+1} = required_ws_vars{k}; %#ok<SAGROW>
    end
end
if ~isempty(missing_ws)
    fprintf('  ERROR: The following variables are missing from the workspace:\n');
    for k = 1:length(missing_ws)
        fprintf('    - %s\n', missing_ws{k});
    end
    fprintf('\n  Re-run run_lbfgsb_optimization_pda_coa.m to regenerate it.\n\n');
    return;
end

% Unpack workspace
opt_outputs = ws.opt_outputs;
clinical    = ws.clinical;
J_final     = ws.J_final;

% Load baseline outputs (pre-opt) for PreOptValue column
if isfield(ws, 'baseline_outputs')
    baseline_outputs = ws.baseline_outputs;
else
    baseline_outputs = [];
end

% Load objective mode
if isfield(ws, 'objective_mode')
    objective_mode = ws.objective_mode;
else
    objective_mode = 'direct_targets_MAP_SV_CoA_mild_zone_P03';
end

% Load CoA penalty settings
coa_penalty_mode    = 'mild_zone';
coa_mild_upper_mmHg = 10.0;
if isfield(ws, 'coa_penalty_mode'),    coa_penalty_mode    = ws.coa_penalty_mode;    end
if isfield(ws, 'coa_mild_upper_mmHg'), coa_mild_upper_mmHg = ws.coa_mild_upper_mmHg; end

% Load objective breakdown for penalty contribution reporting
obj_breakdown_final = [];
if isfield(ws, 'objective_breakdown_final')
    obj_breakdown_final = ws.objective_breakdown_final;
end

if isempty(opt_outputs)
    fprintf('  ERROR: opt_outputs is empty in the workspace.\n');
    fprintf('  The post-optimisation simulation may have failed.\n');
    fprintf('  Re-run run_lbfgsb_optimization_pda_coa.m and check for ODE errors.\n\n');
    return;
end

fprintf('  Patient:         %s\n', safe_get_str(clinical, 'patient_id', 'Unknown'));
fprintf('  Objective mode:  %s\n', objective_mode);
fprintf('  J_final:         %.6f\n\n', J_final);

% =========================================================================
%  STEP 2 — LOAD PATIENT INDEX FROM CSV (for output filename)
% =========================================================================
fprintf('STEP 2: Reading patient data CSV...\n');

csv_path    = fullfile('config', 'patient_data.csv');
patient_idx = NaN;

if ~exist(csv_path, 'file')
    fprintf('  WARNING: patient_data.csv not found. Using patient_id from workspace.\n\n');
    patient_id_str = safe_get_str(clinical, 'patient_id', 'unknown');
else
    try
        patient_table = readtable(csv_path);
        pat_id = safe_get_str(clinical, 'patient_id', '');
        match  = strcmp(patient_table.PatientID, pat_id);
        if any(match)
            patient_idx = find(match, 1);
            fprintf('  Patient %s found at row %d in patient_data.csv.\n\n', pat_id, patient_idx);
        else
            fprintf('  WARNING: Patient ID "%s" not found in CSV.\n\n', pat_id);
        end
        patient_id_str = pat_id;
    catch ME
        fprintf('  WARNING: Could not read patient_data.csv: %s\n\n', ME.message);
        patient_id_str = safe_get_str(clinical, 'patient_id', 'unknown');
    end
end

if ~isnan(patient_idx)
    file_token = sprintf('%d', patient_idx);
else
    file_token = regexprep(patient_id_str, '[^A-Za-z0-9_]', '_');
end

% =========================================================================
%  STEP 3 — ASSEMBLE METRIC DEFINITIONS (direct targets only)
% =========================================================================
fprintf('STEP 3: Assembling metric definitions (MAP, SV, dP_CoA_peak only)...\n');

%  Each row: { MetricLabel, ClinicalField, SimField_pre, SimField_post, Unit }
%
%  Only direct optimization targets are included.
%  SimField_pre  — field name in baseline_outputs (pre-opt simulation)
%  SimField_post — field name in opt_outputs     (post-opt simulation)

metric_defs = {
%  Label          ClinicalField        SimFieldPre        SimFieldPost          Unit
'MAP',          'P_ao_mean_mmHg',    'P_ao_mean',       'P_ao_mean',          'mmHg';
'SV',           'SV_mL',             'SV_lv',           'SV_lv',              'mL';
'dP_CoA_peak',  'dP_coa_mmHg',       'DeltaP_coa_peak', 'DeltaP_coa_peak',    'mmHg';
};

n_metric_defs = size(metric_defs, 1);

% =========================================================================
%  STEP 4 — EXTRACT VALUES AND COMPUTE ERRORS
% =========================================================================
fprintf('STEP 4: Extracting simulated and clinical values...\n\n');

tbl_Metric            = cell(n_metric_defs, 1);
tbl_ClinicalTarget    = NaN(n_metric_defs, 1);
tbl_PreOptValue       = NaN(n_metric_defs, 1);
tbl_PostOptValue      = NaN(n_metric_defs, 1);
tbl_SignedError       = NaN(n_metric_defs, 1);
tbl_AbsError          = NaN(n_metric_defs, 1);
tbl_PercentError      = NaN(n_metric_defs, 1);
tbl_Unit              = cell(n_metric_defs, 1);
tbl_IncludedInObj     = true(n_metric_defs, 1);   % All three are in the objective
tbl_ObjPenaltyMode    = cell(n_metric_defs, 1);   % Penalty mode used in objective
tbl_ObjPenaltyContrib = NaN(n_metric_defs, 1);    % Weighted contribution to J
tbl_Reason            = cell(n_metric_defs, 1);    % Explanation of penalty applied

for k = 1:n_metric_defs
    label         = metric_defs{k, 1};
    clin_field    = metric_defs{k, 2};
    sim_field_pre = metric_defs{k, 3};
    sim_field_pst = metric_defs{k, 4};
    unit_str      = metric_defs{k, 5};

    tbl_Metric{k}  = label;
    tbl_Unit{k}    = unit_str;

    % Clinical value
    clin_val = safe_get_num(clinical, clin_field, NaN);
    tbl_ClinicalTarget(k) = clin_val;

    % Pre-opt simulated value
    pre_val = safe_get_num(baseline_outputs, sim_field_pre, NaN);
    tbl_PreOptValue(k) = pre_val;

    % Post-opt simulated value
    post_val = safe_get_num(opt_outputs, sim_field_pst, NaN);
    tbl_PostOptValue(k) = post_val;

    % Errors: post-opt vs clinical
    if ~isnan(clin_val) && ~isnan(post_val)
        signed_err = post_val - clin_val;
        abs_err    = abs(signed_err);
        if abs(clin_val) > 1e-12
            pct_err = 100 * abs_err / abs(clin_val);
        else
            pct_err = NaN;
        end
        tbl_SignedError(k)  = signed_err;
        tbl_AbsError(k)     = abs_err;
        tbl_PercentError(k) = pct_err;
    end

    % Extract objective penalty breakdown for this metric
    tbl_ObjPenaltyMode{k}    = 'N/A';
    tbl_ObjPenaltyContrib(k) = NaN;
    tbl_Reason{k}            = '';
    if isstruct(obj_breakdown_final) && ~isempty(obj_breakdown_final) && ...
            isfield(obj_breakdown_final, label)
        bd = obj_breakdown_final.(label);
        if isfield(bd, 'weighted_contribution') && isnumeric(bd.weighted_contribution)
            tbl_ObjPenaltyContrib(k) = bd.weighted_contribution;
        end
        if isfield(bd, 'penalty_mode') && (ischar(bd.penalty_mode) || isstring(bd.penalty_mode))
            tbl_ObjPenaltyMode{k} = char(bd.penalty_mode);
        else
            tbl_ObjPenaltyMode{k} = 'standard_normalized';
        end
        if isfield(bd, 'reason') && (ischar(bd.reason) || isstring(bd.reason))
            tbl_Reason{k} = char(bd.reason);
        else
            tbl_Reason{k} = sprintf('Standard normalized squared error (%s)', label);
        end
    end
end

% =========================================================================
%  STEP 5 — COMPUTE SUMMARY STATISTICS (MAP, SV, dP_CoA_peak only)
% =========================================================================
fprintf('STEP 5: Computing summary statistics...\n\n');

% Only include pairs where both clinical and post-opt values are valid
valid_mask = ~isnan(tbl_ClinicalTarget) & ~isnan(tbl_PostOptValue);
n_valid    = sum(valid_mask);

abs_errs_valid = tbl_AbsError(valid_mask);
pct_errs_valid = tbl_PercentError(valid_mask);

if n_valid > 0
    MAE_val              = mean(abs_errs_valid,         'omitnan');
    RMSE_val             = sqrt(mean(abs_errs_valid.^2, 'omitnan'));
    MeanPercentError_val = mean(pct_errs_valid,         'omitnan');
else
    MAE_val              = NaN;
    RMSE_val             = NaN;
    MeanPercentError_val = NaN;
end

% Objective penalty summary (actual J contributions, not clinical error)
obj_valid_mask       = ~isnan(tbl_ObjPenaltyContrib);
n_obj_valid          = sum(obj_valid_mask);
obj_contribs_valid   = tbl_ObjPenaltyContrib(obj_valid_mask);
if n_obj_valid > 0
    TotalObjContrib_val = sum(obj_contribs_valid,  'omitnan');
    MeanObjContrib_val  = mean(obj_contribs_valid, 'omitnan');
else
    TotalObjContrib_val = NaN;
    MeanObjContrib_val  = NaN;
end

% =========================================================================
%  STEP 6 — PRINT CONSOLE TABLE
% =========================================================================
fprintf('=================================================================\n');
fprintf('   DIRECT OPTIMIZATION OUTPUT METRICS — Patient %s\n', ...
    safe_get_str(clinical, 'patient_id', file_token));
fprintf('   Objective mode: %s\n', objective_mode);
fprintf('   CoA penalty:    %s  (mild zone <= %.1f mmHg)\n', coa_penalty_mode, coa_mild_upper_mmHg);
fprintf('   NOTE: For dP_CoA_peak, clinical error is reported vs 4.9 mmHg,\n');
fprintf('   but the objective penalty is zero when sim value <= %.1f mmHg.\n', coa_mild_upper_mmHg);
fprintf('=================================================================\n');

col_w = [14, 12, 12, 12, 11, 11, 11, 7];
hdr = sprintf('%-*s  %*s  %*s  %*s  %*s  %*s  %*s  %*s', ...
    col_w(1), 'Metric', ...
    col_w(2), 'Clinical', ...
    col_w(3), 'PreOpt', ...
    col_w(4), 'PostOpt', ...
    col_w(5), 'SignedErr', ...
    col_w(6), 'AbsErr', ...
    col_w(7), 'PctErr(%)', ...
    col_w(8), 'Unit');
fprintf('%s\n', hdr);
fprintf('%s\n', repmat('-', 1, length(hdr)));

for k = 1:n_metric_defs
    clin_str  = fmt_num(tbl_ClinicalTarget(k), '%10.3f');
    pre_str   = fmt_num(tbl_PreOptValue(k),    '%10.3f');
    post_str  = fmt_num(tbl_PostOptValue(k),   '%10.3f');
    se_str    = fmt_num(tbl_SignedError(k),     '%10.3f');
    ae_str    = fmt_num(tbl_AbsError(k),        '%10.3f');
    pe_str    = fmt_num(tbl_PercentError(k),    '%10.2f');

    fprintf('%-*s  %*s  %*s  %*s  %*s  %*s  %*s  %*s\n', ...
        col_w(1), tbl_Metric{k}, ...
        col_w(2), clin_str, ...
        col_w(3), pre_str, ...
        col_w(4), post_str, ...
        col_w(5), se_str, ...
        col_w(6), ae_str, ...
        col_w(7), pe_str, ...
        col_w(8), tbl_Unit{k});
end

fprintf('%s\n', repmat('=', 1, length(hdr)));

fprintf('\n  CLINICAL AGREEMENT SUMMARY (raw error vs clinical targets):\n');
fprintf('  %-28s  %10s\n', 'Statistic', 'Value');
fprintf('  %s\n', repmat('-', 1, 42));
fprintf('  %-28s  %10.4f\n', 'MAE (mixed units)',  MAE_val);
fprintf('  %-28s  %10.4f\n', 'RMSE (mixed units)', RMSE_val);
fprintf('  %-28s  %10.2f%%\n','Mean Percent Error', MeanPercentError_val);
fprintf('  %-28s  %10d\n',   'N Metrics',          n_valid);
fprintf('  %s\n', repmat('-', 1, 42));

fprintf('\n  OBJECTIVE PENALTY SUMMARY (actual J contributions from optimization):\n');
fprintf('  %-28s  %10s\n', 'Statistic', 'Value');
fprintf('  %s\n', repmat('-', 1, 42));
fprintf('  %-28s  %10.4f\n', 'Total J contribution',   TotalObjContrib_val);
fprintf('  %-28s  %10.4f\n', 'Mean J contribution',    MeanObjContrib_val);
fprintf('  %-28s  %10d\n',   'N active terms',         n_obj_valid);
fprintf('  %s\n', repmat('-', 1, 42));
fprintf('  dP_CoA_peak penalty mode: %s\n', coa_penalty_mode);
fprintf('    Zero penalty if sim dP_CoA_peak <= %.1f mmHg.\n', coa_mild_upper_mmHg);
fprintf('    Clinical error still reported vs %.1f mmHg.\n', ...
    safe_get_num(clinical, 'dP_coa_mmHg', 4.9));

% Show per-metric breakdown
fprintf('\n  Per-metric objective contributions:\n');
for k = 1:n_metric_defs
    if ~isnan(tbl_ObjPenaltyContrib(k))
        fprintf('    %-14s  J_contrib = %.4f  [%s]\n', ...
            tbl_Metric{k}, tbl_ObjPenaltyContrib(k), tbl_Reason{k});
    else
        fprintf('    %-14s  J_contrib = N/A   [%s]\n', tbl_Metric{k}, tbl_ObjPenaltyMode{k});
    end
end

fprintf('\n  Calibration objective J_final = %.6f\n\n', J_final);

% =========================================================================
%  STEP 7 — BUILD AND SAVE OUTPUT TABLES
% =========================================================================
fprintf('STEP 7: Saving CSV tables...\n');

tables_dir = fullfile('results', 'tables');
if ~exist(tables_dir, 'dir')
    mkdir(tables_dir);
    fprintf('  Created directory: %s\n', tables_dir);
end

% --- 7a. Per-metric table ---
T_metrics = table( ...
    tbl_Metric, ...
    tbl_ClinicalTarget, ...
    tbl_PreOptValue, ...
    tbl_PostOptValue, ...
    tbl_SignedError, ...
    tbl_AbsError, ...
    tbl_PercentError, ...
    tbl_Unit, ...
    tbl_IncludedInObj, ...
    tbl_ObjPenaltyMode, ...
    tbl_ObjPenaltyContrib, ...
    tbl_Reason, ...
    'VariableNames', { ...
        'Metric', ...
        'ClinicalTarget', ...
        'PreOptValue', ...
        'PostOptValue', ...
        'SignedError', ...
        'AbsError', ...
        'PercentError', ...
        'Unit', ...
        'IncludedInObjective', ...
        'ObjectivePenaltyMode', ...
        'ObjectivePenaltyContribution', ...
        'Reason' ...
    });

metrics_csv = fullfile(tables_dir, ...
    sprintf('final_validation_metrics_patient_%s.csv', file_token));
writetable(T_metrics, metrics_csv);
fprintf('  Saved per-metric table:  %s\n', metrics_csv);

% --- 7b. Summary table ---
patient_id_cell = {safe_get_str(clinical, 'patient_id', file_token)};
summary_notes   = sprintf(['Direct targets only (MAP,SV,dP_CoA_peak); objective_mode=%s; ' ...
    'coa_penalty=%s(<=%.0fmmHg); J_final=%.6f; N=%d; ' ...
    'For dP_CoA_peak: clinical error vs 4.9mmHg still reported but objective penalty=0 if sim<=%.0fmmHg'], ...
    objective_mode, coa_penalty_mode, coa_mild_upper_mmHg, J_final, n_valid, coa_mild_upper_mmHg);

T_summary = table( ...
    patient_id_cell, ...
    patient_idx, ...
    n_valid, ...
    MAE_val, ...
    RMSE_val, ...
    MeanPercentError_val, ...
    {summary_notes}, ...
    'VariableNames', { ...
        'PatientID', ...
        'PatientIdx', ...
        'NMetrics', ...
        'MAE', ...
        'RMSE', ...
        'MeanPercentError', ...
        'Notes' ...
    });

summary_csv = fullfile(tables_dir, ...
    sprintf('final_validation_summary_patient_%s.csv', file_token));
writetable(T_summary, summary_csv);
fprintf('  Saved summary table:     %s\n', summary_csv);

fprintf('\n');
fprintf('=================================================================\n');
fprintf('   VALIDATION METRICS COMPLETE — Patient %s\n', ...
    safe_get_str(clinical, 'patient_id', file_token));
fprintf('   Metrics: MAP, SV, dP_CoA_peak  |  N = %d\n', n_valid);
fprintf('   Clinical agreement: MAE=%.4f | RMSE=%.4f | MeanPctErr=%.2f%%\n', ...
    MAE_val, RMSE_val, MeanPercentError_val);
fprintf('   Objective penalty:  TotalJ=%.4f | MeanJ=%.4f | CoAPenalty=%s\n', ...
    TotalObjContrib_val, MeanObjContrib_val, coa_penalty_mode);
fprintf('=================================================================\n\n');


%% =========================================================================
%  LOCAL HELPER FUNCTIONS
%% =========================================================================

function val = safe_get_num(s, fname, default_val)
    if isstruct(s) && ~isempty(s) && isfield(s, fname)
        v = s.(fname);
        if isnumeric(v) && isscalar(v)
            val = v;
            return;
        end
    end
    val = default_val;
end


function val = safe_get_str(s, fname, default_val)
    if isstruct(s) && ~isempty(s) && isfield(s, fname)
        v = s.(fname);
        if ischar(v) || isstring(v)
            val = char(v);
            return;
        end
    end
    val = default_val;
end


function str = fmt_num(val, fmt)
    if isnan(val)
        str = 'N/A';
    else
        str = strtrim(sprintf(fmt, val));
    end
end
