# PDA Lumped Parameter Model — Virtual CoA Severity Framework
## Based on: Ortiz-Rangel et al. (2022), Biomedical Signal Processing and Control 71:103151

---

## What This Model Does

This codebase implements a patient-specific Windkessel-type Lumped Parameter Model (LPM) of 
the neonatal cardiovascular system with a Patent Ductus Arteriosus (PDA). It is calibrated 
from PDA patient clinical data and then used to simulate virtual Coarctation of the Aorta 
(CoA) at programmable stenosis severities (50%, 75%, 90% narrowing).

This allows hemodynamic impact assessment of CoA — including pressure gradients, flow 
redistribution, and left ventricular workload — using only PDA-derived physiological data 
as the calibration input. No real CoA patient data is required.

---

## How to Run

1. Ensure `patient_data.csv` is in `config/`
2. Open MATLAB and set the working directory to `pda_lpm/`
3. Run: `main_pda_lpm`
4. Select a patient when prompted
5. Select CoA severity scenarios when prompted

---

## Project Structure

```
pda_lpm/
├── main_pda_lpm.m               ← Entry point. No physics. No plotting.
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
├── tests/
│   ├── test_baseline.m          ← Must pass before any patient run
│   └── test_valve_logic.m       ← Unit test for valve switching
├── results/
│   ├── figures/                 ← PDF exports only
│   └── tables/                  ← CSV summary outputs
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
- CoA extension: P_ao_dist, Q_shunt_coa (12th, 13th states)

---

## References

[1] Ortiz-Rangel E et al. (2022). Dynamic modeling and simulation of the human 
    cardiovascular system with PDA. Biomed Signal Process Control 71:103151.
    https://doi.org/10.1016/j.bspc.2021.103151

[2] Keshavarz-Motamed Z et al. (2011). Effect of coarctation of the aorta and bicuspid 
    aortic valve on flow dynamics and turbulence. J Biomech 44:2817–2825.

[3] Stergiopulos N et al. (1996). Determinants of stroke volume and systolic and 
    diastolic aortic pressure. Am J Physiol 270(6):H2050–H2059.

---

## Author

Cardiovascular Simulation Team  
Date: 2025  
Version: 1.0
