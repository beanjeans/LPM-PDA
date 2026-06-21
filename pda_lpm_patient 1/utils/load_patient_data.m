function [clinical, patient_table] = load_patient_data(csv_path)
% LOAD_PATIENT_DATA
% -----------------------------------------------------------------------
% Loads neonatal patient records from the clinical CSV and returns a
% validated clinical struct for a user-selected patient. All unit
% conversions from raw clinical records are performed here and nowhere else.
%
% INPUTS:
%   csv_path  - path to patient_data.csv                          [string]
%
% OUTPUTS:
%   clinical      - struct of validated clinical measurements for
%                   the selected patient, with fields as specified
%                   in docs/clinical_data_dictionary.md
%   patient_table - full table (all patients) for reference
%
% ASSUMPTIONS:
%   - CSV columns match clinical_data_dictionary.md exactly
%   - Age is in days; TB = Tinggi badan (height) in cm; BB = Berat badan (weight) in grams; BSA in m²
%   - All pressures already in mmHg; velocities in m/s
%
% REFERENCES:
%   [1] docs/clinical_data_dictionary.md
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

%% Load raw table
try
    patient_table = readtable(csv_path);
catch ME
    error('LOAD_PATIENT_DATA: Cannot read %s\n  %s', csv_path, ME.message);
end

n_patients = height(patient_table);

%% Display patient menu
fprintf('\n--- AVAILABLE PATIENTS ---\n');
for k = 1:n_patients
    fprintf('  %d.  %s  |  Age: %3.0f days  |  HR: %3d bpm  |  PDA: %.2f mm\n', ...
        k, ...
        patient_table.PatientID{k}, ...
        patient_table.Age(k), ...
        patient_table.HeartRate(k), ...
        patient_table.DPDA(k));
end
fprintf('\n');

patient_idx = input(sprintf('Select patient number (1-%d): ', n_patients));
if isempty(patient_idx) || patient_idx < 1 || patient_idx > n_patients
    error('LOAD_PATIENT_DATA: Invalid patient number selected.');
end

row = patient_table(patient_idx, :);

%% Build clinical struct — units explicitly converted at load time
% See docs/clinical_data_dictionary.md for field-by-field description

clinical.patient_id         = row.PatientID{1};           % string
clinical.age_days           = row.Age(1);                 % [days]
clinical.sex                = row.Sex{1};                 % 'F' or 'M'

% Body measurements
clinical.height_cm          = row.TB(1);                  % [cm]  — Tinggi badan (body height)
clinical.weight_g           = row.BB(1);                  % [g]   — Berat badan (body weight in grams)
clinical.BSA_m2             = row.BSA(1);                 % [m²]  — directly recorded clinical measurement

% Cardiac timing
clinical.HR_bpm             = row.HeartRate(1);           % [bpm]
clinical.SV_mL              = row.StrokeVolume(1);        % [mL]

% Systemic pressures
clinical.P_ao_sys_mmHg      = row.SSAP(1);               % [mmHg] — systolic
clinical.P_ao_dia_mmHg      = row.SDAP(1);               % [mmHg] — diastolic
clinical.P_ao_mean_mmHg     = row.MAP(1);                 % [mmHg] — mean arterial pressure

% Anatomical diameters (stored in mm; converted to m in physics functions)
clinical.D_shunt_pda_mm     = row.DPDA(1);               % [mm]
clinical.D_coa_mm           = row.DCoA(1);               % [mm] — for reference only
clinical.D_aao_mm           = row.DAAo(1);               % [mm] — ascending aorta
clinical.D_dta_mm           = row.DDTA(1);               % [mm] — descending thoracic aorta
clinical.D_isthmus_mm       = row.DIsthmus(1);           % [mm] — aortic isthmus
clinical.D_dao_mm           = row.DDAo(1);               % [mm] — descending aorta
clinical.D_aov_mm           = row.DAoV(1);               % [mm] — aortic valve annulus
clinical.D_pv_mm            = row.DPV(1);                % [mm] — pulmonary valve annulus

% Doppler velocities
clinical.v_aov_ms           = row.vAoV(1);               % [m/s] — aortic valve (CW Doppler, A5C)
clinical.v_pv_psax_ms       = row.vPVpsax(1);            % [m/s] — pulmonary valve via PSAX (CW Doppler)
clinical.v_pv_supra_ms      = row.vPVsupra(1);           % [m/s] — pulmonary valve via Suprasternal (CW Doppler)
clinical.v_pda_ms           = row.vPDA(1);               % [m/s] — PDA peak velocity (CW Doppler)
clinical.v_coa_ms           = row.vCoA(1);               % [m/s] — CoA peak velocity (ref)

% PDA flow direction (1 = L→R, 2 = R→L, 3 = bidirectional)
clinical.pda_direction      = row.ArahAliranPDA(1);       % [1/2/3]

% Pressure gradients (Doppler-derived, Bernoulli)
clinical.dP_aov_mmHg        = row.dPAoV(1);              % [mmHg]
clinical.dP_pv_mmHg         = row.dPPV(1);               % [mmHg]
clinical.dP_pda_mmHg        = row.dPPDA(1);              % [mmHg]
clinical.dP_coa_mmHg        = row.dPCoA(1);              % [mmHg] — reference for CoA sim

%% Derived quantities (computed once at load time)
% Cardiac output from SV and HR
mLs_to_Lmin             = 60 / 1000;                     % conversion factor [L/min per mL/s]
CO_mLs                  = (clinical.SV_mL * clinical.HR_bpm) / 60; % [mL/s]
clinical.CO_Lmin        = CO_mLs * mLs_to_Lmin;          % [L/min]
clinical.CO_mLs         = CO_mLs;                        % [mL/s]

% Estimated PA pressure from PDA gradient (if L→R shunt: P_ao - P_pa ≈ dP_pda)
if clinical.pda_direction == 1  % Left-to-right
    clinical.P_pa_est_mmHg = max(clinical.P_ao_mean_mmHg - clinical.dP_pda_mmHg, 5);
else
    clinical.P_pa_est_mmHg = clinical.P_ao_mean_mmHg;    % [mmHg] — elevated (R→L or bidirectional)
end

%% Display summary
fprintf('\n--- SELECTED PATIENT: %s ---\n', clinical.patient_id);
fprintf('  Age: %g days  |  Sex: %s  |  Height: %.1f cm  |  Weight: %.0f g\n', ...
    clinical.age_days, clinical.sex, clinical.height_cm, clinical.weight_g);
fprintf('  HR: %d bpm  |  SV: %.2f mL  |  CO: %.2f L/min\n', ...
    clinical.HR_bpm, clinical.SV_mL, clinical.CO_Lmin);
fprintf('  MAP: %.1f mmHg  |  Sys/Dia: %d/%d mmHg\n', ...
    clinical.P_ao_mean_mmHg, clinical.P_ao_sys_mmHg, clinical.P_ao_dia_mmHg);
fprintf('  PDA diameter: %.2f mm  |  PDA direction: %d  |  PDA dP: %.1f mmHg\n', ...
    clinical.D_shunt_pda_mm, clinical.pda_direction, clinical.dP_pda_mmHg);
fprintf('  Estimated PA pressure: %.1f mmHg\n\n', clinical.P_pa_est_mmHg);

end
