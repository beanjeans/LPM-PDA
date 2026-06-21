# PDA Lumped Parameter Model — Virtual CoA Severity Framework
## Based on: Ortiz-Rangel et al. (2022), Biomedical Signal Processing and Control 71:103151

---

## What This Model Does

This codebase implements a patient-specific Windkessel-type Lumped Parameter Model (LPM) of 
the neonatal cardiovascular system with a Patent Ductus Arteriosus (PDA). It is calibrated 
from PDA patient clinical data and then used to simulate virtual Coarctation of the Aorta 
(CoA) at programmable stenosis severities and variable lengths (discrete vs. long-segment).

**Key Features (Version 5.0):**
- **L-BFGS-B Parameter Optimization**: Bounded quasi-Newton calibration of influential parameters against patient clinical targets (MAP, SBP, DBP, SV, CoA gradient), using `fmincon` with physiologically realistic bounds.
- **Global Sensitivity Analysis (GSA)**: A complete Sobol/Saltelli GSA workflow identifies which parameters most influence the CoA clinical outputs, selecting the optimal parameter subset for optimization.
- **Clinical Severity Classification**: Predicts CoA severity (Mild, Moderate, Severe) based on the **simulated pressure gradient**, following ESC guidelines.
- **Variable CoA Length**: Simulates discrete (< 5 mm) and long-segment (≥ 5 mm) coarctations.
- **PDA Confounder Analysis**: Retains PDA parameters to show how a patent PDA can mask observed CoA severity.

This allows hemodynamic impact assessment of CoA — using only PDA-derived physiological data 
as the calibration input. No real CoA patient data is required.

---

## How to Run

> **For the full research pipeline, run the scripts in this exact order:**
> 1. `run_sobol_gsa_pda_coa` → 2. `run_lbfgsb_optimization_pda_coa` → 3. `main_pda_lpm` → 4. `run_final_validation_metrics`

---

### Step 1 — Global Sensitivity Analysis *(run this first)*
Identifies which parameters most influence the CoA outputs.
1. Ensure `patient_data.csv` is in `config/`
2. Open MATLAB and set the working directory to `pda_lpm/`
3. Open `run_sobol_gsa_pda_coa.m` — set `patient_idx` and optionally `N` (default 256).
4. Run: `run_sobol_gsa_pda_coa`
5. Read the final console output — it lists which parameters to include in optimization.
6. Results saved to `results/gsa/`.

### Step 2 — L-BFGS-B Parameter Optimization *(run this second)*
Calibrates the influential parameters to match the patient's clinical measurements.
1. Open `run_lbfgsb_optimization_pda_coa.m` — review **Section A** at the top:
   - Set `patient_idx` to match Step 1.
   - Confirm `opt_param_names` contains the influential parameters from Step 1.
   - Adjust `weights` to prioritize your most trusted clinical targets.
2. Run: `run_lbfgsb_optimization_pda_coa`
3. Note the optimized `stenosis_pct` and `coa_length_mm` printed at the end.
4. Results saved to `results/optimization/`.

### Step 3 — Final Simulation & Figures
Runs the interactive simulation using the calibrated parameter values.
1. Run: `main_pda_lpm`
2. Select the patient when prompted.
3. Enter the **optimized** `stenosis_pct` and `coa_length_mm` from Step 2 when prompted.
4. Generates publication-ready waveform figures and the final CoA severity table.

### Step 4 — Post-Calibration Agreement Metrics *(optional but recommended)*
Reports how well the calibrated model reproduces the clinical training targets.

> **Note**: These are agreement metrics — not independent external validation. The same clinical targets (MAP, SBP, DBP, SV, CoA gradient) were used during calibration.

1. Run: `run_final_validation_metrics`
   - Reads `results/optimization/optimization_workspace.mat` (saved by Step 2).
   - Does **not** re-run GSA, optimization, or simulation — standalone read-only.
2. Outputs saved to `results/tables/`:
   - `final_validation_metrics_patient_<ID>.csv` — per-metric table (SignedError, AbsError, %Error)
   - `final_validation_summary_patient_<ID>.csv` — summary stats (MAE, RMSE, MeanPctErr)
3. Key metrics reported: SBP, DBP, MAP, SV, CoA peak gradient, CO, Qp/Qs.

---

> **`main_pda_lpm` can also be run standalone** (without Steps 1–2) if you simply want
> to explore pressure waveforms at a manually chosen stenosis percentage.

---

## Project Structure

```
pda_lpm/
├── main_pda_lpm.m                     ← Entry point. No physics. No plotting.
├── run_sobol_gsa_pda_coa.m            ← Orchestrator for Global Sensitivity Analysis (GSA).
├── run_lbfgsb_optimization_pda_coa.m  ← Orchestrator for bounded parameter optimization.
├── run_final_validation_metrics.m     ← Post-calibration agreement metrics (Step 4).
├── config/
│   ├── default_parameters.m           ← Reference (healthy neonate) parameter set
│   └── patient_data.csv               ← Clinical records (copy here before running)
├── models/
│   ├── system_rhs_pda.m               ← ODE right-hand side — PDA cardiovascular physics
│   ├── system_rhs_pda_coa.m           ← ODE right-hand side — PDA + virtual CoA physics (13 states)
│   ├── elastance_model.m              ← Piecewise quadratic/cosine elastance (Stergiopulos 1996)
│   └── valve_model.m                  ← Ideal diode valve flow: Q = max(0, dP) / R
├── solvers/
│   └── integrate_system.m             ← Wraps ode15s with documented tolerances
├── utils/
│   ├── unit_conversion.m              ← All conversion factors in one place
│   ├── load_patient_data.m            ← Loads & validates CSV, builds clinical struct
│   │                                     Units: TB [cm], BB [grams], vPV split into
│   │                                     vPVpsax [m/s] and vPVsupra [m/s]
│   ├── build_patient_params.m         ← Calibrates LPM params from clinical data
│   ├── build_coa_params.m             ← Adds virtual CoA geometry + severity thresholds
│   ├── compute_clinical_indices.m     ← Derives: CO, SV, EF, Qp/Qs, SW_lv,
│   │                                     DeltaP_CoA (peak/mean-sys), predicted_CoA_severity
│   └── plot_results.m                 ← Publication-ready figure generation
├── gsa/
│   ├── sample_sobol_params.m          ← Generates Saltelli sampling matrices (8 parameters)
│   ├── evaluate_model_outputs.m       ← Batch model evaluator for GSA
│   └── compute_sobol_indices.m        ← Computes S1/ST indices + 95% bootstrap CI
├── optimization/
│   ├── objective_lbfgsb_pda_coa.m     ← Weighted least-squares objective for fmincon
│   ├── apply_optimized_params.m       ← Converts x_opt vector to params structs
│   └── plot_optimization_results.m    ← Convergence + before/after comparison plots
├── tests/
│   ├── test_baseline.m                ← Must pass before any patient run
│   └── test_valve_logic.m             ← Unit test for valve switching
├── results/
│   ├── figures/                       ← PDF exports only
│   ├── tables/                        ← CSV summary outputs (incl. validation metrics)
│   ├── gsa/                           ← GSA output CSVs, PNGs, and MAT files
│   └── optimization/                  ← Optimization CSVs, PNGs, and MAT workspace
└── docs/
    ├── theory_notes.md                ← Governing equations & assumptions
    └── clinical_data_dictionary.md    ← Maps CSV fields → MATLAB variables
```

---

## Model Equations

The model is a **13-state Windkessel-type LPM** extending Ortiz-Rangel et al. (2022).
See `docs/theory_notes.md` for the full derivation. Plain-text formula reference: [`lpm_formulas.md`](docs/lpm_formulas.md).

**State variables:**

| Group | States |
|-------|--------|
| Base circuit (Ortiz-Rangel 2022) | P_ra, P_rv, P_pa, Q_pa_pul, P_pv, P_la, P_lv, P_ao, Q_ao_sys, P_sys |
| PDA extension | Q_shunt_pda (state 11) |
| CoA extension | P_ao_dist, Q_coa (states 12–13) |

**Key formula building blocks:**

| Component | Formula |
|-----------|---------|
| Valve flow | `Q = max(0, P_up - P_down) / R_valve` |
| Ventricular pressure | `dP_lv/dt = (P_lv/E_lv)*dE_lv/dt + E_lv*(Q_mv - Q_av)` |
| CoA pressure drop | `dP_CoA = R_viscous*Q_CoA + K_turb*Q_CoA*\|Q_CoA\|` |
| Poiseuille resistance | `R_viscous = 128*mu*L / (pi * D^4)` |
| Bernoulli coefficient | `K_turb = 0.5*rho / A_CoA^2` |
| Elastance (contraction) | `En = 1.55 * (tn/Ts1)^2` |
| Elastance (relaxation) | `En = 1.55 * cos(pi/2 * (tn-Ts1)/(Ts2-Ts1))^2` |

**Clinical outputs from `compute_clinical_indices`:**
- Pressures: P_ao sys/dia/mean, P_lv, P_pa, P_ra
- Flows: SV, CO (L/min), Qp/Qs, Q_shunt_PDA
- LV mechanics: EF, stroke work (mmHg·mL and J)
- CoA metrics: DeltaP_CoA peak, mean-systolic gradient, Q_CoA/Q_total
- **Primary output**: `predicted_CoA_severity` — Mild / Moderate / Severe

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
Version: 5.0
