function [params_opt, stenosis_opt, coa_length_opt] = apply_optimized_params(x_opt, opt_config)
% APPLY_OPTIMIZED_PARAMS
% -----------------------------------------------------------------------
% Takes the raw optimizer output vector x_opt and constructs the final
% parameter structs for the post-optimization simulation.
%
% This function is called after fmincon converges to apply the final
% x_opt to the baseline parameter struct and call build_coa_params with
% the optimized stenosis and length values.
%
% DESIGN RULE: This function does NOT modify any existing model file.
% It only assembles the already-computed parameter values into the correct
% structs so that the existing simulate pipeline can be called unchanged.
%
% INPUTS:
%   x_opt       - (D_opt × 1) optimized parameter vector from fmincon
%   opt_config  - same config struct used in objective_lbfgsb_pda_coa
%
% OUTPUTS:
%   params_opt      - complete LPM params struct with optimized values
%                     applied (ready to pass to system_rhs_pda_coa)
%   stenosis_opt    - optimized stenosis_pct value  [%]
%   coa_length_opt  - optimized coa_length_mm value [mm]
%
% AUTHOR:   Optimization Extension — Cardiovascular Simulation Team
% DATE:     2025-01-01
% VERSION:  1.0
% -----------------------------------------------------------------------

%% 1. Start from the fixed baseline params
params_opt     = opt_config.fixed_params;
stenosis_opt   = opt_config.stenosis_pct;    % Default, may be overridden
coa_length_opt = opt_config.coa_length_mm;  % Default, may be overridden

%% 2. Apply each optimized parameter
for k = 1:length(opt_config.param_names)
    pname = opt_config.param_names{k};
    val   = x_opt(k);

    switch pname
        case 'R_systemic'
            params_opt.R_systemic = val;
        case 'C_sys'
            params_opt.C_sys = val;
        case 'Emax_lv'
            params_opt.Emax_lv = val;
            if ~ismember('Emin_lv', opt_config.param_names)
                params_opt.Emin_lv = val * 0.05;
            end
            params_opt.Emax_rv = val * 0.5;
            params_opt.Emin_rv = params_opt.Emin_lv;
        case 'Emin_lv'
            params_opt.Emin_lv = val;
            params_opt.Emin_rv = val;
        case 'stenosis_pct'
            stenosis_opt = val;
        case 'coa_length_mm'
            coa_length_opt = val;
        case 'R_shunt_pda'
            params_opt.R_shunt_pda = val;
        case 'C_ao'
            params_opt.C_ao = val;
        case 'R_pa'
            params_opt.R_pa = val;
    end
end

%% 3. Clamp to hard physical limits
% (should already be within bounds from fmincon, but guard defensively)
params_opt.R_systemic  = max(0.01, params_opt.R_systemic);
params_opt.C_sys       = max(1e-4, params_opt.C_sys);
params_opt.Emax_lv     = max(0.01, params_opt.Emax_lv);
params_opt.Emin_lv     = max(1e-4, params_opt.Emin_lv);
params_opt.Emax_rv     = max(0.01, params_opt.Emax_rv);
params_opt.Emin_rv     = max(1e-4, params_opt.Emin_rv);
params_opt.C_ao        = max(1e-4, params_opt.C_ao);
params_opt.R_pa        = max(0.001, params_opt.R_pa);
params_opt.R_shunt_pda = max(0.01, params_opt.R_shunt_pda);
stenosis_opt           = max(0.1, min(99.9, stenosis_opt));
coa_length_opt         = max(0.1, coa_length_opt);

fprintf('  [apply_optimized_params] Parameters applied successfully.\n');

end
