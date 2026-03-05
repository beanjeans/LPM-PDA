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
DeltaP_coa = R_coa_viscous * Q_coa + K_turb * Q_coa * |Q_coa|

dP_ao_dist/dt = (Q_coa - Q_ao_sys_dist) / C_sys_dist
dQ_coa/dt = (P_ao - P_ao_dist - DeltaP_coa) / L_coa
```

---

## 2. Elastance Model (Stergiopulos 1996)

The normalized elastance En(tn) follows a piecewise quadratic / cosine model 
adapted from the original Stergiopulos formulation used in Ortiz-Rangel (2022):

- Contraction phase (0 ≤ tn ≤ Ts1):  En = 1.55 * (tn/Ts1)²
- Relaxation phase (Ts1 < tn ≤ Ts2): En = 1.55 * cos²(π/2 * (tn-Ts1)/(Ts2-Ts1))
- Diastole (tn > Ts2):                En = 0

where Ts1 = 0.3*sqrt(T), Ts2 = 0.45*sqrt(T), T = cardiac cycle period.

---

## 3. Valve Logic

Heart valves modelled as ideal diodes with resistance:

  Q_valve = max(0, P_upstream - P_downstream) / R_valve

The `max(0, ...)` form avoids hard switching discontinuities (guardrail §8.4).

---

## 4. CoA Stenosis Geometry

Stenosis percentage S is defined as:
  S = (1 - (D_stenosis / D_reference)²) × 100%

where D_stenosis is effective orifice diameter and D_reference is normal aortic diameter.

Given S:
  D_stenosis_m = D_reference_m * sqrt(1 - S/100)
  A_coa_m2     = π * (D_stenosis_m/2)²

Viscous resistance (Poiseuille):
  R_coa_viscous = (128 * μ * L_coa) / (π * D_stenosis_m⁴)   [Pa·s/m³]
  → convert to mmHg·s/mL via unit_conversion factors

Turbulent coefficient (Bernoulli):
  K_turb = 0.5 * ρ / A_coa_m²   [Pa/(m³/s)²]
  → convert to mmHg/(mL/s)²

---

## 5. Assumptions

1. Blood is incompressible and Newtonian (viscosity μ = 0.004 Pa·s)
2. Lumped-parameter (0D) — spatial gradients neglected
3. PDA shunt is purely left-to-right at baseline (validated by P_ao > P_pa for neonates)
4. CoA is modelled as a fixed geometric obstruction (no wall compliance at CoA)
5. All chamber compliances are constant (no series compliance model)
6. Neonate scaling applied to all compliances and initial conditions
