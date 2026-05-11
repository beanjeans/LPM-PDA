function uc = unit_conversion()
% UNIT_CONVERSION
% -----------------------------------------------------------------------
% Returns a struct of all unit conversion factors used in the PDA LPM.
% Centralising conversions eliminates silent unit mixing across files.
%
% OUTPUTS:
%   uc  - struct of named, commented conversion factors [dimensionless ratios]
%
% USAGE:
%   uc = unit_conversion();
%   P_Pa = P_mmHg * uc.mmHg_to_Pa;
%
% REFERENCES:
%   [1] NIST SI unit definitions (physics.nist.gov)
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

% --- Pressure ---
uc.mmHg_to_Pa  = 133.322;    % [Pa / mmHg]  — 1 mmHg = 133.322 Pa (NIST)
uc.Pa_to_mmHg  = 1 / 133.322;

% --- Volume / Flow ---
uc.m3s_to_mLs  = 1e6;        % [mL/s per m³/s]
uc.mLs_to_m3s  = 1e-6;       % [m³/s per mL/s]
uc.mLs_to_Lmin = 60 / 1000;  % [L/min per mL/s]
uc.Lmin_to_mLs = 1000 / 60;  % [mL/s per L/min]

% --- Length / Area ---
uc.mm_to_m     = 1e-3;        % [m / mm]
uc.m_to_mm     = 1e3;         % [mm / m]
uc.cm_to_m     = 1e-2;        % [m / cm]
uc.m_to_cm     = 1e2;         % [cm / m]
uc.mm2_to_m2   = 1e-6;        % [m² / mm²]
uc.m2_to_mm2   = 1e6;         % [mm² / m²]

% --- Resistance unit conversion ---
% Internal model units: mmHg·s/mL
% Fluid-mechanics formula (Poiseuille) produces: Pa·s/m³
%   R_mmHg_s_mL = R_Pa_s_m3 * Pa_to_mmHg * mLs_to_m3s
uc.Pa_s_m3_to_mmHg_s_mL = uc.Pa_to_mmHg * uc.mLs_to_m3s;  % ≈ 7.5e-9

% --- Turbulent coefficient unit conversion ---
% K_turb formula produces: Pa / (m³/s)²
%   K_mmHg_mLs2 = K_Pa_m3s2 * Pa_to_mmHg * (mLs_to_m3s)^2
uc.Pa_m3s2_to_mmHg_mLs2 = uc.Pa_to_mmHg * (uc.mLs_to_m3s)^2;  % ≈ 7.5e-21

% --- Inertance unit conversion ---
% L formula produces: Pa·s²/m³
%   L_mmHg_s2_mL = L_Pa_s2_m3 * Pa_to_mmHg * mLs_to_m3s
uc.Pa_s2_m3_to_mmHg_s2_mL = uc.Pa_to_mmHg * uc.mLs_to_m3s;  % same as resistance factor

% --- Energy ---
uc.mmHg_mL_to_J = 1e-3 / 7.5006;  % [J / (mmHg·mL)] — 1 mmHg·mL = 1.333e-4 J
uc.J_to_mmHg_mL = 1 / uc.mmHg_mL_to_J;

% --- Time ---
uc.min_to_s    = 60;           % [s / min]
uc.s_to_min    = 1 / 60;

end
