# Theory Notes — PDA LPM with Virtual CoA Framework

## 1. Governing Equations

### 1.1 Normal Cardiovascular Circuit (Ortiz-Rangel 2022, Eqs. 8–17)

State variables and their ODEs (pressure in mmHg, flow in mL/s, volume in mL):

```
dP_ra/dt = (P_sys - P_ra)/(R_systemic * C_ra) + D_tv*(P_rv - P_ra)/(R_tv * C_ra)
dP_rv/dt = D_tv*(P_ra - P_rv)/(R_tv*C_rv) + D_pv*(P_pa - P_rv)/(R_pv*C_rv) - P_rv*dE_rv/E_rv
dP_pa/dt = D_pv*(P_rv - P_pa)/(R_pv*C_pa) - Q_pa_pul/C_pa
dQ_pa_pul/dt = (P_pa - P_pv)/L_pa - Q_pa_pul*R_pa/L_pa
dP_pv/dt = Q_pa_pul/C_pv + (P_la - P_pv)/(R_pv_veins*C_pv)
dP_la/dt = (P_pv - P_la)/(R_pv_veins*C_la) + D_mv*(P_lv - P_la)/(R_mv*C_la)
dP_lv/dt = D_mv*(P_la - P_lv)/(R_mv*C_lv) + D_av*(P_ao - P_lv)/(R_av*C_lv) - P_lv*dE_lv/E_lv
dP_ao/dt = D_av*(P_lv - P_ao)/(R_av*C_ao) - Q_ao_sys/C_ao
dQ_ao_sys/dt = (P_ao - P_sys)/L_ao - Q_ao_sys*R_ao/L_ao
dP_sys/dt = Q_ao_sys/C_sys + (P_ra - P_sys)/(R_systemic*C_sys)
```

### 1.2 PDA Extension (Ortiz-Rangel 2022, Eqs. 18–20)

The PDA connects the aorta to the pulmonary artery with resistance R_shunt_pda and
inertance L_shunt_pda:

```
dP_pa/dt += Q_shunt_pda / C_pa          (additional inflow to PA)
dP_ao/dt += Q_shunt_pda / C_ao          (additional outflow from Ao)
dQ_shunt_pda/dt = (P_pa - P_ao)/L_shunt_pda - Q_shunt_pda*R_shunt_pda/L_shunt_pda
```

Sign convention: Q_shunt_pda > 0 = left-to-right (aorta → PA), physiological direction.

### 1.3 Virtual CoA Extension

CoA is modelled as an added resistance R_coa_viscous (Poiseuille) plus a quadratic
turbulent term K_turb in the aortic segment, splitting the aorta into proximal (P_ao)
and distal (P_ao_dist) compartments:

```
ΔP_coa = R_coa_viscous × Q_coa + K_turb × Q_coa × |Q_coa|

dP_ao_dist/dt = (Q_coa - Q_sys_lower) / C_sys
dQ_coa/dt     = (P_ao - P_ao_dist - ΔP_coa) / L_coa
```

Both R_coa_viscous and K_turb (Bernoulli coefficient) depend on the effective CoA
orifice area derived from `stenosis_pct`, and on the CoA segment length `coa_length_mm`.
See Section 4 for geometry details.

---

## 2. Elastance Model (Stergiopulos 1996)

The normalized elastance En(tn) follows a piecewise quadratic / cosine model
adapted from the original Stergiopulos formulation used in Ortiz-Rangel (2022):

- Contraction phase (0 ≤ tn ≤ Ts1):  En = 1.55 × (tn/Ts1)²
- Relaxation phase (Ts1 < tn ≤ Ts2): En = 1.55 × cos²(π/2 × (tn−Ts1)/(Ts2−Ts1))
- Diastole (tn > Ts2):                En = 0

where Ts1 = 0.3×√T, Ts2 = 0.45×√T, T = cardiac cycle period.

---

## 3. Valve Logic

Heart valves modelled as ideal diodes with resistance:

  Q_valve = max(0, P_upstream − P_downstream) / R_valve

The `max(0, ...)` form avoids hard switching discontinuities (guardrail §8.4).

---

## 4. CoA Stenosis Geometry

Stenosis percentage S is defined as area reduction:

  S = (1 − (A_stenosis / A_reference)) × 100%
    = (1 − (D_stenosis / D_reference)²) × 100%

where D_stenosis is effective orifice diameter and D_reference is normal aortic diameter.

Given S and D_reference:
```
D_stenosis_m = D_reference_m × sqrt(1 − S/100)
A_coa_m²     = π × (D_stenosis_m/2)²
```

Viscous resistance (Poiseuille), Ref [2] Eq. 3:
```
R_coa_viscous = (128 × μ × L_coa) / (π × D_stenosis_m⁴)   [Pa·s/m³]
              → convert to mmHg·s/mL via unit_conversion factors
```

Turbulent coefficient (Bernoulli), Ref [2] Eq. 4–5:
```
K_turb = 0.5 × ρ / A_coa_m²   [Pa/(m³/s)²]
        → convert to mmHg/(mL/s)²
```

CoA segment inertance, Ref [2]:
```
L_coa = (4 × ρ × L_coa) / (π × D_stenosis_m²)   [Pa·s²/m³]
       → convert to mmHg·s²/mL
```

**Note**: Both R_coa_viscous and L_coa scale with `coa_length_mm`; longer CoA segments
produce higher viscous resistance and inertance, amplifying the pressure gradient.

---

## 5. CoA Length Classification

CoA segment length should **not** be assumed as a fixed constant.
Patient-specific length from echocardiography, CT, or MRI should be used when available.

When patient data are unavailable, simulate across the following canonical scenarios:

| Scenario             | Length      | Category              |
|----------------------|-------------|-----------------------|
| Discrete CoA         | 3 mm        | discrete/short        |
| Threshold case       | 5 mm        | boundary              |
| Long-segment CoA     | 8–10 mm     | long-segment          |

**Literature classification threshold** (Ref [4]):
- **Discrete / short CoA**: length < 5 mm
- **Long-segment CoA**: length ≥ 5 mm

This boundary is defined in `build_coa_params.m` as `COA_LENGTH_THRESHOLD_MM = 5`.
Edit that variable there to change the category boundary.

---

## 6. CoA Severity Classification (Primary Clinical Output)

CoA severity is **not** reported directly from `stenosis_pct`.
The anatomical stenosis drives CoA geometry and the simulated pressure drop;
clinical severity is then classified from that simulated pressure gradient,
consistent with clinical practice (Doppler/cath gradient classification).

### 6.1 Gradient Metrics

For each simulation, three gradient metrics are extracted from the last
steady-state cardiac cycle:

| Metric                   | Symbol              | Description                                   |
|--------------------------|---------------------|-----------------------------------------------|
| Peak gradient            | ΔP_CoA_peak         | max(P_ao − P_ao_dist) over the cycle          |
| Mean systolic gradient   | ΔP_CoA_mean_sys     | mean gradient during LV ejection phase        |
| Full-cycle mean          | ΔP_CoA_mean         | mean(P_ao − P_ao_dist) over entire cycle      |

**Classification** uses `ΔP_CoA_mean_sys` (mean-systolic gradient), following Ref [3]:

| Classification   | Mean Systolic ΔP_CoA   |
|------------------|------------------------|
| **Mild**         | < 20 mmHg              |
| **Moderate**     | 20–40 mmHg             |
| **Severe**       | > 40 mmHg              |

These thresholds are stored in `params.severity_thresholds` (set in `build_coa_params`)
and are editable without changing any other file:

```matlab
sev_thr.mild_upper_mmHg     = 20;   % edit here only
sev_thr.moderate_upper_mmHg = 40;   % edit here only
params_coa.severity_thresholds = sev_thr;
```

### 6.2 Output Fields per Simulation

For each CoA scenario, `compute_clinical_indices` returns:

| Field                    | Units    | Description                                     |
|--------------------------|----------|-------------------------------------------------|
| `stenosis_pct`           | %        | Anatomical input (geometry driver only)         |
| `coa_length_mm`          | mm       | CoA segment length used                         |
| `coa_length_category`    | string   | 'discrete/short' or 'long-segment'              |
| `DeltaP_coa_peak`        | mmHg     | Peak pressure gradient across CoA               |
| `DeltaP_coa_mean_sys`    | mmHg     | Mean systolic gradient (used for classification)|
| `DeltaP_coa_mean`        | mmHg     | Full-cycle mean gradient                        |
| `Q_coa_fraction`         | [0–1]    | Q_coa_mean / Q_total_mean                       |
| `P_ao_proximal_mean`     | mmHg     | Mean proximal aortic pressure                   |
| `P_ao_distal_mean`       | mmHg     | Mean distal aortic pressure                     |
| `predicted_CoA_severity` | string   | 'mild', 'moderate', or 'severe'                 |
| `pda_modifier_note`      | string   | Note on PDA masking effect                      |

---

## 7. PDA as Haemodynamic Modifier

A patent ductus arteriosus (PDA) can **mask or modify** the pressure gradient
across CoA and must be reported as a **confounder**, not as part of CoA anatomy:

1. **Retrograde PDA flow** (right-to-left or bidirectional):  
   Blood from PA enters the descending aorta distal to the CoA, raising P_ao_dist.  
   This reduces ΔP_CoA → may cause **under-classification** of CoA severity.

2. **Left-to-right PDA shunt**:  
   Aortic blood diverts to PA, reducing flow through CoA → lowers Q_coa.  
   This may paradoxically reduce ΔP_CoA despite high anatomical stenosis.

**Implementation**: PDA parameters (`R_shunt_pda`, `L_shunt_pda`) are always retained
in `params_coa`. The field `model.pda_modifier_note` in the output reports whether a
PDA is present and warns that the observed ΔP_CoA may underestimate true severity.

---

## 8. Assumptions

1. Blood is incompressible and Newtonian (viscosity μ = 0.004 Pa·s)
2. Lumped-parameter (0D) — spatial gradients neglected
3. PDA shunt is primarily left-to-right at baseline (P_ao > P_pa for neonates)
4. CoA is modelled as a fixed geometric obstruction (no wall compliance at CoA)
5. All chamber compliances are constant (no series compliance model)
6. Neonate scaling applied to all compliances and initial conditions
7. CoA length is **not** fixed; it must be supplied as `coa_length_mm` in
   `build_coa_params`. Use scenario values (3 / 5 / 8–10 mm) when patient
   measurement is unavailable.
8. Predicted severity is based on the **simulated mean-systolic ΔP_CoA**, not
   on `stenosis_pct` directly, to reflect clinical diagnostic practice.

---

## 9. References

| #   | Citation |
|-----|----------|
| [1] | Ortiz-Rangel et al. (2022). *Biomed Signal Process Control* 71:103151. |
| [2] | Keshavarz-Motamed et al. (2011). *J Biomech* 44:2817–2825. |
| [3] | Baumgartner H et al. (2010). *Eur Heart J* 31(19):2369–2417. (ESC valvular/CoA guidelines: significant gradient ≥ 20 mmHg; severe > 40 mmHg) |
| [4] | Vergales JE et al. (2013). *Pediatr Cardiol* 34:1616–1623. / Campbell M et al. (2002). (CoA length classification: discrete < 5 mm; long-segment ≥ 5 mm) |
| [5] | Stergiopulos N et al. (1996). *Am J Physiol* 270(6):H2050–H2059. (Elastance normalisation) |
| [6] | Rudolph AM (2001). *Congenital Diseases of the Heart*. Futura Pub. (Neonate haemodynamic reference values) |
