# Clinical Data Dictionary

Maps fields in `patient_data.csv` to MATLAB variables used in the LPM.

| CSV Field     | MATLAB Variable              | Unit (CSV) | Unit (Model) | Method            | Reliability | Notes                              |
|---------------|------------------------------|------------|--------------|-------------------|-------------|------------------------------------|
| PatientID     | clinical.patient_id          | string     | string       | Record            | High        | De-identified                      |
| Age           | clinical.age_days            | days       | days         | Record            | High        |                                    |
| Sex           | clinical.sex                 | F/M        | F/M          | Record            | High        |                                    |
| TB            | clinical.weight_total_g      | g          | g            | Scale             | High        | Total body weight                  |
| BB            | clinical.weight_body_kg      | kg         | kg           | Scale             | High        | Body weight (non-fluid)            |
| BSA           | clinical.BSA_m2              | m²         | m²           | Clinical record   | Moderate    | Directly recorded; not formula-derived |
| HeartRate     | clinical.HR_bpm              | bpm        | bpm          | ECG               | High        |                                    |
| StrokeVolume  | clinical.SV_mL               | mL         | mL           | Echo              | Moderate    |                                    |
| BP            | clinical.BP_string           | mmHg       | —            | Sphygmomanometer  | Moderate    | Formatted as "sys/dia"             |
| SSAP          | clinical.P_ao_sys_mmHg       | mmHg       | mmHg         | Sphygmomanometer  | Moderate    | Systolic aortic pressure           |
| SDAP          | clinical.P_ao_dia_mmHg       | mmHg       | mmHg         | Sphygmomanometer  | Moderate    | Diastolic aortic pressure          |
| MAP           | clinical.P_ao_mean_mmHg      | mmHg       | mmHg         | Derived           | Derived     | Mean arterial pressure             |
| DPDA          | clinical.D_shunt_pda_mm      | mm         | mm → m (internal) | Echo        | Moderate    | PDA duct diameter                  |
| DCoA          | clinical.D_coa_mm            | mm         | mm → m (internal) | Echo        | Moderate    | CoA orifice diameter (if present)  |
| DAAo          | clinical.D_aao_mm            | mm         | mm           | Echo              | Moderate    | Ascending aorta diameter           |
| DDTA          | clinical.D_dta_mm            | mm         | mm           | Echo              | Moderate    | Descending thoracic aorta          |
| DIsthmus      | clinical.D_isthmus_mm        | mm         | mm           | Echo              | Moderate    | Aortic isthmus diameter            |
| DDAo          | clinical.D_dao_mm            | mm         | mm           | Echo              | Moderate    | Descending aorta diameter          |
| DAoV          | clinical.D_aov_mm            | mm         | mm           | Echo              | Moderate    | Aortic valve annulus               |
| DPV           | clinical.D_pv_mm             | mm         | mm           | Echo              | Moderate    | Pulmonary valve annulus            |
| vAoV          | clinical.v_aov_ms            | m/s        | m/s          | Doppler           | High        | Peak velocity through aortic valve |
| vPV           | clinical.v_pv_ms             | m/s        | m/s          | Doppler           | High        | Peak velocity through pulm valve   |
| vPDA          | clinical.v_pda_ms            | m/s        | m/s          | Doppler           | High        | Peak PDA flow velocity             |
| vCoA          | clinical.v_coa_ms            | m/s        | m/s          | Doppler           | High        | Peak CoA velocity (if present)     |
| ArahAliranPDA | clinical.pda_direction       | 1/2/3      | 1/2/3        | Doppler           | High        | 1=L→R, 2=R→L, 3=bidirectional     |
| dPAoV         | clinical.dP_aov_mmHg         | mmHg       | mmHg         | Doppler (Bernoulli)| Moderate   | Gradient across aortic valve       |
| dPPV          | clinical.dP_pv_mmHg          | mmHg       | mmHg         | Doppler (Bernoulli)| Moderate   | Gradient across pulm valve         |
| dPPDA         | clinical.dP_pda_mmHg         | mmHg       | mmHg         | Doppler (Bernoulli)| Moderate   | PDA pressure gradient              |
| dPCoA         | clinical.dP_coa_mmHg         | mmHg       | mmHg         | Doppler (Bernoulli)| Moderate   | CoA gradient (reference only)      |

## Notes on Reliability

- **High**: Direct measurement, small inter-observer variability
- **Moderate**: Echo-derived, ±10–20% typical uncertainty
- **Derived**: Computed from other measurements, error propagates from source measurements

## Unit Conversions Applied at Load Time

All unit conversions are performed in `utils/load_patient_data.m`:
- Diameters: mm → stored as mm; converted to m only inside physics functions
- BSA: already in m² in the CSV (directly recorded clinical measurement; no formula applied)
- EF: not directly in CSV; derived in `compute_clinical_indices.m`
