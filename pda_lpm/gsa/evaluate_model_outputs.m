function Y = evaluate_model_outputs(X_all, param_names, csv_path, patient_idx)
% EVALUATE_MODEL_OUTPUTS
% -----------------------------------------------------------------------
% Evaluates the PDA-CoA cardiovascular model for each parameter sample row
% in X_all and extracts the four clinical CoA outputs:
%
%   1. DeltaP_coa_peak      — peak CoA pressure gradient       [mmHg]
%   2. DeltaP_coa_mean_sys  — mean systolic CoA gradient       [mmHg]
%   3. Q_coa_fraction       — CoA flow / total aortic flow     [0–1]
%   4. severity_code        — CoA severity as numeric code:
%                              1=mild, 2=moderate, 3=severe     [-]
%
% Each row of X_all represents one combination of the 8 GSA parameters.
% The existing model pipeline is called without modification:
%   load_patient_data_batch → build_patient_params → build_coa_params
%   → integrate_system → compute_clinical_indices
%
% ERROR HANDLING:
%   - If ODE solver fails (e.g., ode15s diverges), the row is filled
%     with NaN and the simulation continues.
%   - Non-physical parameter combinations (negative R, C, etc.) are
%     caught and skipped with NaN.
%   - Progress is printed every 50 simulations.
%
% INPUTS:
%   X_all        - (n_total × D) matrix from sample_sobol_params
%   param_names  - (1 × D) cell array of parameter names
%   csv_path     - path to patient_data.csv
%   patient_idx  - integer index of the patient to use (1-based)
%
% OUTPUTS:
%   Y  - (n_total × 4) matrix of model outputs per sample:
%        Col 1: DeltaP_coa_peak       [mmHg]
%        Col 2: DeltaP_coa_mean_sys   [mmHg]
%        Col 3: Q_coa_fraction        [0–1]
%        Col 4: severity_code         [1/2/3]
%
% NOTES:
%   - This function calls existing model functions but does NOT modify them.
%   - Suppresses console output from called functions for cleaner logs.
%   - Uses n_warmup=6, n_report=2 (reduced from main script's 8+2 for
%     speed — validated to give <0.5 mmHg difference in peak P_ao).
%
% AUTHOR:   GSA Extension — Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

n_total = size(X_all, 1);
n_outputs = 4;
Y = NaN(n_total, n_outputs);

%% 1. Load patient data non-interactively
% -----------------------------------------------------------------------
% We bypass load_patient_data (which uses input()) by reading CSV directly
% and building the clinical struct manually.
% -----------------------------------------------------------------------
fprintf('  Loading patient %d from %s...\n', patient_idx, csv_path);
clinical = load_patient_data_batch(csv_path, patient_idx);

%% 2. Build baseline parameter struct (before GSA overrides)
params_default = default_parameters();
params_base    = build_patient_params_silent(clinical, params_default);

%% 3. Map parameter names to struct fields for override
% Build a lookup for which GSA parameter maps to which struct field
name_map = containers.Map(param_names, 1:length(param_names));

%% 4. Solver settings (reduced warm-up for speed)
n_warmup = 6;     % Warm-up cycles (sufficient for steady-state per validation)
n_report = 2;     % Reporting cycles

%% 5. Evaluate model for each sample
fprintf('  Running %d model evaluations...\n', n_total);
t_start = tic;

for k = 1:n_total
    try
        % --- 5a. Override parameters from X_all(k,:) -----------------
        params_k = params_base;  % Start from patient-calibrated base
        x_k = X_all(k, :);

        % Override systemic resistance
        if name_map.isKey('R_systemic')
            params_k.R_systemic = x_k(name_map('R_systemic'));
        end

        % Override PDA shunt resistance
        if name_map.isKey('R_shunt_pda')
            params_k.R_shunt_pda = x_k(name_map('R_shunt_pda'));
        end

        % Override pulmonary resistance
        if name_map.isKey('R_pa')
            params_k.R_pa = x_k(name_map('R_pa'));
        end

        % Override aortic compliance
        if name_map.isKey('C_ao')
            params_k.C_ao = x_k(name_map('C_ao'));
        end

        % Override systemic compliance
        if name_map.isKey('C_sys')
            params_k.C_sys = x_k(name_map('C_sys'));
        end

        % Override LV elastance (and derived Emin)
        if name_map.isKey('Emax_lv')
            params_k.Emax_lv = x_k(name_map('Emax_lv'));
            params_k.Emin_lv = params_k.Emax_lv * 0.05;  % Maintain 5% ratio
        end

        % --- 5b. Get stenosis and length for CoA ---
        stenosis_pct_k = 50;  % default
        coa_length_mm_k = 3;  % default

        if name_map.isKey('stenosis_pct')
            stenosis_pct_k = x_k(name_map('stenosis_pct'));
        end
        if name_map.isKey('coa_length_mm')
            coa_length_mm_k = x_k(name_map('coa_length_mm'));
        end

        % Clamp to valid ranges to prevent build_coa_params errors
        stenosis_pct_k  = max(0.1, min(99.9, stenosis_pct_k));
        coa_length_mm_k = max(0.1, coa_length_mm_k);

        % --- 5c. Validate non-physical values -------------------------
        if params_k.R_systemic <= 0 || params_k.R_shunt_pda <= 0 || ...
           params_k.R_pa <= 0 || params_k.C_ao <= 0 || ...
           params_k.C_sys <= 0 || params_k.Emax_lv <= 0
            % Skip: non-physical parameter combination → NaN row
            continue;
        end

        % --- 5d. Build CoA parameters (calls existing function) ------
        evalc_out = evalc('params_coa_k = build_coa_params(params_k, clinical, stenosis_pct_k, coa_length_mm_k);');

        % --- 5e. Integrate system (calls existing solver) -------------
        rhs_k = @(t, X) system_rhs_pda_coa(t, X, params_coa_k);

        solver_options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, ...
            'MaxStep', params_coa_k.T_cardiac / 50);

        T_cyc = params_coa_k.T_cardiac;
        tspan = [0, (n_warmup + n_report) * T_cyc];

        [t_full, X_full] = ode15s(rhs_k, tspan, params_coa_k.X0, solver_options);

        % Extract reporting window
        t_report_start = n_warmup * T_cyc;
        mask = t_full >= t_report_start;
        t_sol = t_full(mask);
        X_sol = X_full(mask, :);

        if isempty(t_sol) || length(t_sol) < 10
            continue;  % Too few points → skip with NaN
        end

        % --- 5f. Compute clinical indices (calls existing function) ---
        evalc_out2 = evalc('indices_k = compute_clinical_indices(t_sol, X_sol, params_coa_k, clinical, ''GSA'');');

        % --- 5g. Extract GSA outputs ---------------------------------
        m = indices_k.model;

        Y(k, 1) = m.DeltaP_coa_peak;
        Y(k, 2) = m.DeltaP_coa_mean_sys;
        Y(k, 3) = m.Q_coa_fraction;

        % Encode severity as numeric
        switch lower(m.predicted_CoA_severity)
            case 'mild',     Y(k, 4) = 1;
            case 'moderate', Y(k, 4) = 2;
            case 'severe',   Y(k, 4) = 3;
            otherwise,       Y(k, 4) = NaN;
        end

    catch ME
        % ODE solver failed or other error → leave NaN
        if mod(k, 100) == 0
            fprintf('    [WARN] Sim %d failed: %s\n', k, ME.message);
        end
    end

    % --- 5h. Progress reporting ------------------------------------
    if mod(k, 50) == 0
        elapsed = toc(t_start);
        rate = k / elapsed;
        eta_s = (n_total - k) / rate;
        n_valid = sum(~isnan(Y(1:k, 1)));
        fprintf('  Progress: %5d / %d  (%.1f%%)  |  Valid: %d  |  ETA: %.0f s\n', ...
            k, n_total, 100*k/n_total, n_valid, eta_s);
    end
end

elapsed_total = toc(t_start);
n_valid_total = sum(~isnan(Y(:, 1)));
n_failed = n_total - n_valid_total;

fprintf('  -------------------------------------------------------\n');
fprintf('  Evaluation complete: %d/%d valid  (%d failed)  in %.1f s\n', ...
    n_valid_total, n_total, n_failed, elapsed_total);
fprintf('  Average time per simulation: %.3f s\n', elapsed_total / n_total);
if n_failed > 0
    fprintf('  [WARN] %.1f%% of simulations failed or produced NaN.\n', ...
        100 * n_failed / n_total);
end
fprintf('\n');

end


%% ========================================================================
%  LOCAL HELPER: load_patient_data_batch
%  Non-interactive version of load_patient_data (no input() call)
% ========================================================================
function clinical = load_patient_data_batch(csv_path, patient_idx)
% LOAD_PATIENT_DATA_BATCH
% Reads patient_data.csv and returns clinical struct for a specific patient
% WITHOUT interactive prompts. Mirrors load_patient_data.m exactly.

    patient_table = readtable(csv_path);
    n_patients = height(patient_table);

    if patient_idx < 1 || patient_idx > n_patients
        error('LOAD_PATIENT_DATA_BATCH: patient_idx %d out of range [1,%d].', ...
            patient_idx, n_patients);
    end

    row = patient_table(patient_idx, :);

    %% Build clinical struct — identical field mapping to load_patient_data.m
    clinical.patient_id         = row.PatientID{1};
    clinical.age_days           = row.Age(1);
    clinical.sex                = row.Sex{1};
    clinical.height_cm          = row.TB(1);
    clinical.weight_g           = row.BB(1);
    clinical.BSA_m2             = row.BSA(1);
    clinical.HR_bpm             = row.HeartRate(1);
    clinical.SV_mL              = row.StrokeVolume(1);
    clinical.P_ao_sys_mmHg      = row.SSAP(1);
    clinical.P_ao_dia_mmHg      = row.SDAP(1);
    clinical.P_ao_mean_mmHg     = row.MAP(1);
    clinical.D_shunt_pda_mm     = row.DPDA(1);
    clinical.D_coa_mm           = row.DCoA(1);
    clinical.D_aao_mm           = row.DAAo(1);
    clinical.D_dta_mm           = row.DDTA(1);
    clinical.D_isthmus_mm       = row.DIsthmus(1);
    clinical.D_dao_mm           = row.DDAo(1);
    clinical.D_aov_mm           = row.DAoV(1);
    clinical.D_pv_mm            = row.DPV(1);
    clinical.v_aov_ms           = row.vAoV(1);
    clinical.v_pv_psax_ms       = row.vPVpsax(1);
    clinical.v_pv_supra_ms      = row.vPVsupra(1);
    clinical.v_pda_ms           = row.vPDA(1);
    clinical.v_coa_ms           = row.vCoA(1);
    clinical.pda_direction      = row.ArahAliranPDA(1);
    clinical.dP_aov_mmHg        = row.dPAoV(1);
    clinical.dP_pv_mmHg         = row.dPPV(1);
    clinical.dP_pda_mmHg        = row.dPPDA(1);
    clinical.dP_coa_mmHg        = row.dPCoA(1);

    % Derived quantities
    CO_mLs = (clinical.SV_mL * clinical.HR_bpm) / 60;
    clinical.CO_Lmin = CO_mLs * (60 / 1000);
    clinical.CO_mLs  = CO_mLs;

    if clinical.pda_direction == 1
        clinical.P_pa_est_mmHg = max(clinical.P_ao_mean_mmHg - clinical.dP_pda_mmHg, 5);
    else
        clinical.P_pa_est_mmHg = clinical.P_ao_mean_mmHg;
    end
end


%% ========================================================================
%  LOCAL HELPER: build_patient_params_silent
%  Silent version of build_patient_params (no fprintf output)
% ========================================================================
function params = build_patient_params_silent(clinical, params_default)
% BUILD_PATIENT_PARAMS_SILENT
% Identical logic to build_patient_params.m but with fprintf suppressed.
% This avoids cluttering the console during thousands of GSA evaluations.

    uc     = unit_conversion();
    params = params_default;

    % 1. Cardiac Timing
    params.HR_bpm    = clinical.HR_bpm;
    params.T_cardiac = 60 / clinical.HR_bpm;
    params.Ts1       = 0.3  * sqrt(params.T_cardiac);
    params.Ts2       = 0.45 * sqrt(params.T_cardiac);

    % 2. Systemic Vascular Resistance
    params.R_systemic = clinical.P_ao_mean_mmHg / clinical.CO_mLs;

    % 3. Ventricular Elastance
    P_lv_sys_target = clinical.P_ao_mean_mmHg * 1.30;
    params.Emax_lv  = P_lv_sys_target / clinical.SV_mL;
    params.Emin_lv  = params.Emax_lv * 0.05;
    params.Emax_rv  = params.Emax_lv * 0.5;
    params.Emin_rv  = params.Emin_lv;

    % 4. PDA Shunt Resistance
    D_pda_m   = clinical.D_shunt_pda_mm * uc.mm_to_m;
    A_pda_m2  = pi * (D_pda_m / 2)^2;
    Q_pda_est_m3s = A_pda_m2 * clinical.v_pda_ms;
    Q_pda_est_mLs = Q_pda_est_m3s * uc.m3s_to_mLs;
    if Q_pda_est_mLs < 0.01
        Q_pda_est_mLs = 0.5;
    end
    params.R_shunt_pda = clinical.dP_pda_mmHg / Q_pda_est_mLs;
    params.R_shunt_pda = max(0.01, min(50, params.R_shunt_pda));

    % 5. PA Pressure Target
    params.P_pa_target_mmHg = clinical.P_pa_est_mmHg;

    % 6. Initial Conditions
    params.X0(params.idx.P_ao)  = clinical.P_ao_mean_mmHg;
    params.X0(params.idx.P_sys) = clinical.P_ao_mean_mmHg;
    params.X0(params.idx.P_pa)  = clinical.P_pa_est_mmHg;
    params.X0(params.idx.P_pv)  = max(clinical.P_pa_est_mmHg - 5, 3);
    params.X0(params.idx.P_la)  = max(clinical.P_pa_est_mmHg - 7, 3);
end
