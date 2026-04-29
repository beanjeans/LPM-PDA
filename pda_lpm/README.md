# PDA Lumped Parameter Model — Virtual CoA Severity Framework
## Based on: Ortiz-Rangel et al. (2022), Biomedical Signal Processing and Control 71:103151

---

## What This Model Does

This codebase implements a patient-specific Windkessel-type Lumped Parameter Model (LPM) of 
the neonatal cardiovascular system with a Patent Ductus Arteriosus (PDA). It is calibrated 
from PDA patient clinical data and then used to simulate virtual Coarctation of the Aorta 
(CoA) at programmable stenosis severities and variable lengths (discrete vs. long-segment).

**Key Features (Version 3.0):**
- **Global Sensitivity Analysis (GSA)**: Includes a complete Sobol/Saltelli GSA workflow to identify which model parameters most strongly influence clinical CoA outputs, helping select parameters for optimization.
- **Clinical Severity Classification**: Predicts CoA severity (Mild, Moderate, Severe) based on the **simulated pressure gradient**, following ESC guidelines, rather than just raw anatomical stenosis.
- **Variable CoA Length**: Simulates discrete (< 5 mm) and long-segment (≥ 5 mm) coarctations.
- **PDA Confounder Analysis**: Retains PDA parameters to demonstrate how a patent PDA can mask the observed CoA pressure gradient and lead to severity underestimation.

This allows hemodynamic impact assessment of CoA — using only PDA-derived physiological data 
as the calibration input. No real CoA patient data is required.

---

## How to Run

**Baseline Simulation:**
1. Ensure `patient_data.csv` is in `config/`
2. Open MATLAB and set the working directory to `pda_lpm/`
3. Run: `main_pda_lpm`
4. Select a patient when prompted
5. Select CoA stenosis percentages and segment lengths when prompted

**Global Sensitivity Analysis:**
1. Open MATLAB and set the working directory to `pda_lpm/`
2. Open `run_sobol_gsa_pda_coa.m` to configure the base sample size `N` and `patient_idx`.
3. Run: `run_sobol_gsa_pda_coa`
4. The script will generate sampling matrices, evaluate the model in batch, and output results to `results/gsa/`.

---

## Project Structure

```
pda_lpm/
├── main_pda_lpm.m               ← Entry point. No physics. No plotting.
├── run_sobol_gsa_pda_coa.m      ← Orchestrator for Global Sensitivity Analysis (GSA).
├── config/
│   ├── default_parameters.m     ← Reference (healthy neonate) parameter set
│   └── patient_data.csv         ← Clinical records (copy here before running)
├── models/
│   ├── system_rhs_pda.m         ← ODE right-hand side — PDA cardiovascular physics
│   ├── system_rhs_pda_coa.m     ← ODE right-hand side — PDA + virtual CoA physics
│   ├── elastance_model.m        ← Time-varying ventricular/atrial elastance
│   └── valve_model.m            ← Valve and PDA shunt flow logic
├── solvers/
│   └── integrate_system.m       ← Wraps ode15s with documented tolerances
├── utils/
│   ├── unit_conversion.m        ← All conversion factors in one place
│   ├── load_patient_data.m      ← Loads & validates CSV, builds clinical struct
│   ├── build_patient_params.m   ← Calibrates LPM params from clinical data
│   ├── build_coa_params.m       ← Adds virtual CoA parameters to param struct
│   ├── compute_clinical_indices.m ← Derived hemodynamic indices: CO, SV, Qp/Qs, SW
│   └── plot_results.m           ← Publication-ready figure generation
├── gsa/
│   ├── sample_sobol_params.m    ← Generates Saltelli sampling matrices
│   ├── evaluate_model_outputs.m ← Batch model evaluator for GSA
│   └── compute_sobol_indices.m  ← Computes S1/ST indices and confidence intervals
├── tests/
│   ├── test_baseline.m          ← Must pass before any patient run
│   └── test_valve_logic.m       ← Unit test for valve switching
├── results/
│   ├── figures/                 ← PDF exports only
│   ├── tables/                  ← CSV summary outputs
│   └── gsa/                     ← GSA output CSVs, PNGs, and MAT files
└── docs/
    ├── theory_notes.md          ← Governing equations & assumptions
    └── clinical_data_dictionary.md ← Maps CSV fields → MATLAB variables
```

---

## Model Equations

The model follows the 10-state Windkessel LPM of Ortiz-Rangel et al. (2022).
See `docs/theory_notes.md` for full equation derivation.

Core states:
- P_ra, P_rv, P_pa, Q_pa_pul, P_pv, P_la, P_lv, P_ao, Q_ao_sys, P_sys
- PDA extension: Q_shunt_pda (11th state)
- CoA extension: P_ao_dist, Q_coa (12th, 13th states)

---

## References

[1] Ortiz-Rangel E et al. (2022). Dynamic modeling and simulation of the human 
    cardiovascular system with PDA. Biomed Signal Process Control 71:103151.
    https://doi.org/10.1016/j.bspc.2021.103151

[2] Keshavarz-Motamed Z et al. (2011). Effect of coarctation of the aorta and bicuspid 
    aortic valve on flow dynamics and turbulence. J Biomech 44:2817–2825.

[3] Baumgartner H et al. (2010). ESC Guidelines for the management of grown-up 
    congenital heart disease. Eur Heart J 31(19):2369–2417.

[4] Vergales JE et al. (2013). Native neonatal coarctation of the aorta: length 
    predicts intervention. Pediatr Cardiol 34:1616–1623.

[5] Stergiopulos N et al. (1996). Determinants of stroke volume and systolic and 
    diastolic aortic pressure. Am J Physiol 270(6):H2050–H2059.

---

## Author

Author
Cardiovascular Simulation Team  
Date: 2025  
Version: 3.0
