%% RUN_FINAL_VALIDATION_METRICS
% =========================================================================
% FINAL VALIDATION METRICS — PDA-CoA LPM — Patient 1 (PDA-only)
%
% PURPOSE:
%   Compute and report agreement between the final post-optimisation
%   simulation outputs and the patient's clinical reference values,
%   restricted to the direct optimization targets only.
%
%   Patient 1 is PDA-only.  The optimization objective uses only MAP and SV.
%   Final validation metrics are therefore computed only for MAP and SV.
%   CoA gradient, SBP, DBP, CO, and Qp/Qs are NOT included because:
%     - dP_CoA is not measured for this patient (CSV value = 0, not clinical 0)
%     - SBP, DBP are not direct targets in this 3-parameter setup
%     - CO and Qp/Qs were not optimized
%
% WORKFLOW POSITION:
%   Run AFTER run_lbfgsb_optimization_pda_coa.m has been completed.
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
% METRICS REPORTED:
%   MAP, SV  (direct optimization targets only)
%
% SUMMARY STATISTICS (over MAP and SV only):
%   MAE, RMSE, MeanPercentError
%
% RULES:
%   - Does NOT modify any existing .m files.
%   - Uses NaN for missing values; excludes NaN pairs from summary stats.
%   - Division-by-zero in percent error is handled safely (returns NaN).
%
% REFERENCES:
%   [1] Ortiz-Rangel et al. (2022). Biomed Signal Process Control 71:103151.
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  2.0  — restricted to MAP+SV for PDA-only patient
% =========================================================================

clear; clc;

fprintf('=================================================================\n');
fprintf('   FINAL VALIDATION METRICS — DIRECT OPTIMIZATION OUTPUTS\n');
fprintf('   PDA-CoA Lumped Parameter Cardiovascular Model\n');
fprintf('=================================================================\n\n');
fprintf('  Patient 1 is PDA-only.  Optimization objective uses MAP and SV only.\n');
fprintf('  Final validation metrics are computed only for direct optimization\n');
fprintf('  outputs: MAP and SV.\n');
fprintf('  This optimization uses only MAP and SV as direct clinical targets.\n');
fprintf('  Other model outputs are not included in the optimization or final\n');
fprintf('  validation metrics for this configuration.\n\n');

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

% Verify expected variables exist in the workspace
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
    fprintf('\n  The workspace may be from an older pipeline version.\n');
    fprintf('  Re-run run_lbfgsb_optimization_pda_coa.m to regenerate it.\n\n');
    return;
end

% Unpack relevant workspace contents
opt_outputs = ws.opt_outputs;   % Post-optimisation simulated outputs struct
clinical    = ws.clinical;      % Clinical reference struct
J_final     = ws.J_final;       % Final objective value

if isempty(opt_outputs)
    fprintf('  ERROR: opt_outputs is empty in the workspace.\n');
    fprintf('  The post-optimisation simulation may have failed.\n');
    fprintf('  Re-run run_lbfgsb_optimization_pda_coa.m and check for ODE errors.\n\n');
    return;
end

fprintf('  Patient: %s\n', safe_get_str(clinical, 'patient_id', 'Unknown'));
fprintf('  Final optimisation objective J = %.6f\n\n', J_final);

% =========================================================================
%  STEP 2 — LOAD PATIENT INDEX FROM CSV (for output filename)
% =========================================================================
fprintf('STEP 2: Reading patient data CSV...\n');

csv_path = fullfile('config', 'patient_data.csv');
patient_idx = NaN;   % Default if CSV unavailable

if ~exist(csv_path, 'file')
    fprintf('  WARNING: patient_data.csv not found at: %s\n', csv_path);
    fprintf('  Will use patient_id from workspace for output filenames.\n\n');
    patient_id_str = safe_get_str(clinical, 'patient_id', 'unknown');
else
    try
        patient_table = readtable(csv_path);
        % Find which row matches the loaded clinical patient_id
        pat_id = safe_get_str(clinical, 'patient_id', '');
        match  = strcmp(patient_table.PatientID, pat_id);
        if any(match)
            patient_idx = find(match, 1);
            fprintf('  Patient %s found at row %d in patient_data.csv.\n\n', pat_id, patient_idx);
        else
            fprintf('  WARNING: Patient ID "%s" not found in CSV. Using ID string for filename.\n\n', pat_id);
        end
        patient_id_str = pat_id;
    catch ME
        fprintf('  WARNING: Could not read patient_data.csv: %s\n\n', ME.message);
        patient_id_str = safe_get_str(clinical, 'patient_id', 'unknown');
    end
end

% Build a clean filename token from the patient ID / index
if ~isnan(patient_idx)
    file_token = sprintf('%d', patient_idx);
else
    file_token = regexprep(patient_id_str, '[^A-Za-z0-9_]', '_');
end

% =========================================================================
%  STEP 3 — ASSEMBLE METRIC DEFINITIONS
% =========================================================================
fprintf('STEP 3: Assembling metric definitions...\n');

%  Each row: { MetricLabel, ClinicalField, SimField, Unit, InOptimization, Notes }
%
%  Only MAP and SV are included — these are the direct optimization targets
%  for Patient 1 (PDA-only).  SBP, DBP, dP_CoA, CO, Qp/Qs are excluded:
%    - dP_CoA: not measured for Patient 1 (CSV = 0, not a clinical 0 mmHg)
%    - SBP, DBP: not direct targets in this 3-parameter setup
%    - CO, Qp/Qs: not optimization targets

metric_defs = {
%  Label    ClinicalField        SimField       Unit    InOpt    Notes
'MAP',    'P_ao_mean_mmHg',    'P_ao_mean',   'mmHg',  'true',  'Mean arterial pressure; direct target (w=5.0)';
'SV',     'SV_mL',             'SV_lv',       'mL',    'true',  'Stroke volume; direct target (w=4.0)';
};

n_metric_defs = size(metric_defs, 1);

% =========================================================================
%  STEP 4 — EXTRACT VALUES AND COMPUTE ERRORS
% =========================================================================
fprintf('STEP 4: Extracting simulated and clinical values...\n\n');

% Pre-allocate output table arrays
tbl_Metric              = cell(n_metric_defs, 1);
tbl_ClinicalTarget      = NaN(n_metric_defs, 1);
tbl_SimulatedValue      = NaN(n_metric_defs, 1);
tbl_SignedError         = NaN(n_metric_defs, 1);
tbl_AbsoluteError       = NaN(n_metric_defs, 1);
tbl_PercentError        = NaN(n_metric_defs, 1);
tbl_Unit                = cell(n_metric_defs, 1);
tbl_IncludedInOpt       = cell(n_metric_defs, 1);
tbl_Notes               = cell(n_metric_defs, 1);

for k = 1:n_metric_defs
    label         = metric_defs{k, 1};
    clin_field    = metric_defs{k, 2};
    sim_field     = metric_defs{k, 3};
    unit_str      = metric_defs{k, 4};
    in_opt_str    = metric_defs{k, 5};
    notes_str     = metric_defs{k, 6};

    tbl_Metric{k}        = label;
    tbl_Unit{k}          = unit_str;
    tbl_IncludedInOpt{k} = in_opt_str;
    tbl_Notes{k}         = notes_str;

    % --- Retrieve clinical value ---
    if isempty(clin_field)
        clin_val = NaN;   % No clinical target defined for this metric
    else
        clin_val = safe_get_num(clinical, clin_field, NaN);
    end
    tbl_ClinicalTarget(k) = clin_val;

    % --- Retrieve simulated value ---
    sim_val = safe_get_num(opt_outputs, sim_field, NaN);
    tbl_SimulatedValue(k) = sim_val;

    % --- Compute error metrics (only when both values are available) ---
    if ~isnan(clin_val) && ~isnan(sim_val)
        signed_err  = sim_val - clin_val;
        abs_err     = abs(signed_err);
        if abs(clin_val) > 1e-12          % Guard against divide-by-zero
            pct_err = 100 * abs_err / abs(clin_val);
        else
            pct_err = NaN;                % Clinical value is zero — percent error undefined
        end
        tbl_SignedError(k)   = signed_err;
        tbl_AbsoluteError(k) = abs_err;
        tbl_PercentError(k)  = pct_err;
    else
        % Missing pair — errors remain NaN; excluded from summary stats
        tbl_SignedError(k)   = NaN;
        tbl_AbsoluteError(k) = NaN;
        tbl_PercentError(k)  = NaN;
    end
end

% =========================================================================
%  STEP 5 — COMPUTE SUMMARY STATISTICS
% =========================================================================
fprintf('STEP 5: Computing summary statistics...\n\n');

% Only include pairs where BOTH clinical and simulated values are valid
valid_mask  = ~isnan(tbl_ClinicalTarget) & ~isnan(tbl_SimulatedValue);
n_valid     = sum(valid_mask);

abs_errs_valid = tbl_AbsoluteError(valid_mask);
pct_errs_valid = tbl_PercentError(valid_mask);

if n_valid > 0
    MAE_val              = mean(abs_errs_valid,         'omitnan');
    RMSE_val             = sqrt(mean(abs_errs_valid.^2, 'omitnan'));
    MeanPercentError_val = mean(pct_errs_valid,         'omitnan');
    MaxAbsoluteError_val = max(abs_errs_valid,          [], 'omitnan');
else
    MAE_val              = NaN;
    RMSE_val             = NaN;
    MeanPercentError_val = NaN;
    MaxAbsoluteError_val = NaN;
end

% =========================================================================
%  STEP 6 — PRINT CONSOLE TABLE
% =========================================================================
fprintf('=================================================================\n');
fprintf('   FINAL VALIDATION METRICS — Patient %s\n', ...
    safe_get_str(clinical, 'patient_id', file_token));
fprintf('   Direct optimization outputs: MAP and SV only.\n');
fprintf('=================================================================\n');

col_w = [14, 14, 14, 11, 11, 11, 7];   % Column widths
hdr = sprintf('%-*s  %*s  %*s  %*s  %*s  %*s  %*s', ...
    col_w(1), 'Metric', ...
    col_w(2), 'Clinical', ...
    col_w(3), 'Simulated', ...
    col_w(4), 'SignedErr', ...
    col_w(5), 'AbsErr', ...
    col_w(6), 'PctErr(%)', ...
    col_w(7), 'Unit');
fprintf('%s\n', hdr);
fprintf('%s\n', repmat('-', 1, length(hdr)));

for k = 1:n_metric_defs
    clin_str = fmt_num(tbl_ClinicalTarget(k), '%10.3f');
    sim_str  = fmt_num(tbl_SimulatedValue(k), '%10.3f');
    se_str   = fmt_num(tbl_SignedError(k),    '%10.3f');
    ae_str   = fmt_num(tbl_AbsoluteError(k),  '%10.3f');
    pe_str   = fmt_num(tbl_PercentError(k),   '%10.2f');

    fprintf('%-*s  %*s  %*s  %*s  %*s  %*s  %*s\n', ...
        col_w(1), tbl_Metric{k}, ...
        col_w(2), clin_str, ...
        col_w(3), sim_str, ...
        col_w(4), se_str, ...
        col_w(5), ae_str, ...
        col_w(6), pe_str, ...
        col_w(7), tbl_Unit{k});
end

fprintf('%s\n', repmat('=', 1, length(hdr)));
fprintf('\n  Summary Statistics (MAP and SV only — %d direct optimization outputs):\n', n_valid);
fprintf('  %-28s  %10s\n', 'Statistic', 'Value');
fprintf('  %s\n', repmat('-', 1, 42));
fprintf('  %-28s  %10.4f\n', 'MAE (over MAP and SV)',    MAE_val);
fprintf('  %-28s  %10.4f\n', 'RMSE (over MAP and SV)',   RMSE_val);
fprintf('  %-28s  %10.2f%%\n','Mean Percent Error',       MeanPercentError_val);
fprintf('  %-28s  %10d\n',   'N Metrics included',        n_valid);
fprintf('  %s\n', repmat('=', 1, 42));
fprintf('\n  Final optimization objective J_final = %.6f\n', J_final);
fprintf('  No measured CoA target is used for this PDA-only patient.\n\n');

% =========================================================================
%  STEP 7 — BUILD AND SAVE OUTPUT TABLES
% =========================================================================
fprintf('STEP 7: Saving CSV tables...\n');

% --- Create output directory if needed ---
tables_dir = fullfile('results', 'tables');
if ~exist(tables_dir, 'dir')
    mkdir(tables_dir);
    fprintf('  Created directory: %s\n', tables_dir);
end

% --- 7a. Per-metric table ---
T_metrics = table( ...
    tbl_Metric, ...
    tbl_ClinicalTarget, ...
    tbl_SimulatedValue, ...
    tbl_SignedError, ...
    tbl_AbsoluteError, ...
    tbl_PercentError, ...
    tbl_Unit, ...
    tbl_IncludedInOpt, ...
    tbl_Notes, ...
    'VariableNames', { ...
        'Metric', ...
        'ClinicalTarget', ...
        'PreOptValue', ...
        'SignedError', ...
        'AbsError', ...
        'PercentError', ...
        'Unit', ...
        'IncludedInObjective', ...
        'Notes' ...
    });

metrics_csv = fullfile(tables_dir, ...
    sprintf('final_validation_metrics_patient_%s.csv', file_token));
writetable(T_metrics, metrics_csv);
fprintf('  Saved per-metric table:  %s\n', metrics_csv);

% --- 7b. Summary table ---
patient_id_cell  = {safe_get_str(clinical, 'patient_id', file_token)};
summary_notes    = sprintf('Direct optimization outputs MAP+SV only; J_final=%.6f; N_valid=%d', ...
                           J_final, n_valid);

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
fprintf('   FINAL VALIDATION COMPLETE\n');
fprintf('   Patient: %s  |  Direct optimization outputs: MAP, SV\n', ...
    safe_get_str(clinical, 'patient_id', file_token));
fprintf('   MAE = %.4f  |  RMSE = %.4f  |  MeanPctErr = %.2f%%\n', ...
    MAE_val, RMSE_val, MeanPercentError_val);
fprintf('=================================================================\n\n');


%% =========================================================================
%  LOCAL HELPER FUNCTIONS
%% =========================================================================

function val = safe_get_num(s, fname, default_val)
% SAFE_GET_NUM  Safely retrieve a numeric field from a struct.
%   Returns default_val if the struct is empty, the field does not exist,
%   or the field value is not a finite scalar number.
    if isstruct(s) && isfield(s, fname)
        v = s.(fname);
        if isnumeric(v) && isscalar(v)
            val = v;
            return;
        end
    end
    val = default_val;
end


function val = safe_get_str(s, fname, default_val)
% SAFE_GET_STR  Safely retrieve a character/string field from a struct.
%   Returns default_val string if the struct is empty, the field does not
%   exist, or the field value is not a character/string.
    if isstruct(s) && isfield(s, fname)
        v = s.(fname);
        if ischar(v) || isstring(v)
            val = char(v);
            return;
        end
    end
    val = default_val;
end


function str = fmt_num(val, fmt)
% FMT_NUM  Format a number using fmt, or return 'N/A' if NaN.
    if isnan(val)
        str = 'N/A';
    else
        str = sprintf(fmt, val);
        str = strtrim(str);
    end
end
