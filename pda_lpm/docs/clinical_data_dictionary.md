# Clinical Data Dictionary

Maps fields in `patient_data.csv` to MATLAB variables used in the LPM.

## Actively Used Fields

These fields are read by the model and directly used in physics calculations,
parameter calibration, or the optimisation objective function.

| CSV Field     | MATLAB Variable           | Unit (CSV) | Unit (Model)      | Used In                                          |
|---------------|---------------------------|------------|-------------------|--------------------------------------------------|
| PatientID     | clinical.patient_id       | string     | string            | Patient identification, output filenames         |
| HeartRate     | clinical.HR_bpm           | bpm        | bpm               | Cardiac timing (`T_cardiac`, `Ts1`, `Ts2`)       |
| StrokeVolume  | clinical.SV_mL            | mL         | mL                | Elastance calibration (`Emax_lv`), CO, obj fn    |
| SSAP          | clinical.P_ao_sys_mmHg    | mmHg       | mmHg              | Optimisation objective (SBP target)              |
| SDAP          | clinical.P_ao_dia_mmHg    | mmHg       | mmHg              | Optimisation objective (DBP target)              |
| MAP           | clinical.P_ao_mean_mmHg   | mmHg       | mmHg              | `R_systemic`, initial conditions, obj fn (MAP)   |
| DPDA          | clinical.D_shunt_pda_mm   | mm         | mm → m (internal) | `Q_pda_est` → `R_shunt_pda` computation         |
| DAAo          | clinical.D_aao_mm         | mm         | mm → m (internal) | Reference aortic diameter for virtual CoA geometry|
| vPDA          | clinical.v_pda_ms         | m/s        | m/s               | `Q_pda_est` → `R_shunt_pda` computation         |
| ArahAliranPDA | clinical.pda_direction    | 1/2/3      | 1/2/3             | Estimated PA pressure (`P_pa_est_mmHg`)          |
| dPPDA         | clinical.dP_pda_mmHg      | mmHg       | mmHg              | `R_shunt_pda`, PA pressure est., obj fn          |
| dPCoA         | clinical.dP_coa_mmHg      | mmHg       | mmHg              | Optimisation objective (CoA gradient target)     |

## Stored-Only Fields (Reference Data — Not Used in Calculations)

These fields are loaded into the `clinical` struct but are **not referenced**
in any physics function, calibration step, or objective function.
They are retained as clinical reference data only.

| CSV Field  | MATLAB Variable          | Unit   | Notes                                         |
|------------|--------------------------|--------|-----------------------------------------------|
| Age        | clinical.age_days        | days   | Patient descriptor; printed to console only   |
| Sex        | clinical.sex             | F/M    | Patient descriptor; printed to console only   |
| TB         | clinical.height_cm        | cm     | Tinggi badan (body height); patient descriptor; not used in model |
| BB         | clinical.weight_g         | g      | Berat badan (body weight in grams); patient descriptor; not used in model |
| BSA        | clinical.BSA_m2          | m²     | Directly recorded clinical value; not used in scaling |
| BP         | clinical.BP_string       | mmHg   | Raw string ("sys/dia"); superseded by SSAP/SDAP |
| DCoA       | clinical.D_coa_mm        | mm     | Echo CoA orifice; stored for reference only   |
| DDTA       | clinical.D_dta_mm        | mm     | Descending thoracic aorta; not in physics     |
| DIsthmus   | clinical.D_isthmus_mm    | mm     | Aortic isthmus; not in physics                |
| DDAo       | clinical.D_dao_mm        | mm     | Descending aorta; not in physics              |
| DAoV       | clinical.D_aov_mm        | mm     | Aortic valve annulus; not in physics          |
| DPV        | clinical.D_pv_mm         | mm     | Pulmonary valve annulus; not in physics       |
| vAoV       | clinical.v_aov_ms        | m/s    | AoV peak velocity; not in physics             |
| vPVpsax    | clinical.v_pv_psax_ms     | m/s    | PV peak velocity via PSAX (CW Doppler); not in physics         |
| vPVsupra   | clinical.v_pv_supra_ms    | m/s    | PV peak velocity via Suprasternal (CW Doppler); not in physics |
| vCoA       | clinical.v_coa_ms        | m/s    | CoA peak velocity; stored for reference only  |
| dPAoV      | clinical.dP_aov_mmHg     | mmHg   | AoV gradient; not in physics or obj fn        |
| dPPV       | clinical.dP_pv_mmHg      | mmHg   | PV gradient; not in physics or obj fn         |

## Notes on Reliability

- **High**: Direct measurement, small inter-observer variability
- **Moderate**: Echo-derived, ±10–20% typical uncertainty
- **Derived**: Computed from other measurements, error propagates from source measurements

## Unit Conversions Applied at Load Time

All unit conversions are performed in `utils/load_patient_data.m`:
- Diameters: mm → stored as mm; converted to m only inside physics functions
- TB (Tinggi badan): stored as-is in cm; no conversion applied
- BB (Berat badan): stored as-is in grams; no conversion applied
- BSA: already in m² in the CSV (directly recorded clinical measurement; no formula applied)
- EF: not directly in CSV; derived in `compute_clinical_indices.m`
