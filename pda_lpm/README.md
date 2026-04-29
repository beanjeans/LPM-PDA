# PDA Lumped Parameter Model — Virtual CoA Severity Framework
## Based on: Ortiz-Rangel et al. (2022), Biomedical Signal Processing and Control 71:103151

---

## What This Model Does

This codebase implements a patient-specific Windkessel-type Lumped Parameter Model (LPM) of 
the neonatal cardiovascular system with a Patent Ductus Arteriosus (PDA). It is calibrated 
from PDA patient clinical data and then used to simulate virtual Coarctation of the Aorta 
(CoA) at programmable stenosis severities and variable lengths (discrete vs. long-segment).

**Key Features (Version 4.0):**
- **L-BFGS-B Parameter Optimization**: Bounded quasi-Newton calibration of influential parameters against patient clinical targets (MAP, SBP, DBP, SV, CoA gradient), using `fmincon` with physiologically realistic bounds.
- **Global Sensitivity Analysis (GSA)**: A complete Sobol/Saltelli GSA workflow identifies which parameters most influence the CoA clinical outputs, selecting the optimal parameter subset for optimization.
- **Clinical Severity Classification**: Predicts CoA severity (Mild, Moderate, Severe) based on the **simulated pressure gradient**, following ESC guidelines.
- **Variable CoA Length**: Simulates discrete (< 5 mm) and long-segment (≥ 5 mm) coarctations.
- **PDA Confounder Analysis**: Retains PDA parameters to show how a patent PDA can mask observed CoA severity.

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
2. Open `run_sobol_gsa_pda_coa.m` to configure `N` and `patient_idx`.
3. Run: `run_sobol_gsa_pda_coa`
4. Outputs go to `results/gsa/`. The final console report lists influential parameters.

**L-BFGS-B Parameter Optimization:**
1. Open MATLAB and set the working directory to `pda_lpm/`
2. Open `run_lbfgsb_optimization_pda_coa.m` — configure `patient_idx`, `opt_param_names`, `opt_bounds`, and `weights` in Section A.
3. Run: `run_lbfgsb_optimization_pda_coa`
4. Outputs go to `results/optimization/`: `optimized_parameters.csv`, `objective_history.csv`, and three PNGs.

---

## Project Structure

```
pda_lpm/
├── main_pda_lpm.m               ← Entry point. No physics. No plotting.
├── run_sobol_gsa_pda_coa.m      ← Orchestrator for Global Sensitivity Analysis (GSA).
├── run_lbfgsb_optimization_pda_coa.m  ← Orchestrator for bounded parameter optimization.
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
├── optimization/
│   ├── objective_lbfgsb_pda_coa.m  ← Weighted objective function for fmincon
│   ├── apply_optimized_params.m    ← Converts x_opt vector to params structs
│   └── plot_optimization_results.m ← Convergence + before/after comparison plots
├── tests/
│   ├── test_baseline.m          ← Must pass before any patient run
│   └── test_valve_logic.m       ← Unit test for valve switching
├── results/
│   ├── figures/                 ← PDF exports only
│   ├── tables/                  ← CSV summary outputs
│   ├── gsa/                     ← GSA output CSVs, PNGs, and MAT files
│   └── optimization/            ← Optimization CSVs, PNGs, and MAT workspace
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

Cardiovascular Simulation Team  
Date: 2025  
Version: 4.0
