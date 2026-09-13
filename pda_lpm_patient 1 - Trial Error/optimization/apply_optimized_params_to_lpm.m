function [params, applied_list, ignored_list] = apply_optimized_params_to_lpm(params, optimized_table, disease_mode)
% APPLY_OPTIMIZED_PARAMS_TO_LPM
% -----------------------------------------------------------------------
% Applies optimized parameter values from a CSV-loaded table to the LPM
% parameter struct, respecting disease mode restrictions.
%
% Called by main_pda_lpm.m as the fallback when optimization_workspace.mat
% is unavailable but optimized_parameters.csv is present.
%
% For PDA_only patients: applies Emax_lv and R_systemic only.
%   stenosis_pct, coa_length_mm, R_pa, and R_shunt_pda are ignored.
%   Emin_lv, Emax_rv, Emin_rv are derived from the optimized Emax_lv value.
%
% For PDA_CoA patients: applies all available parameters including
%   stenosis_pct, R_pa, R_shunt_pda if present in the table.
%
% INPUTS:
%   params          - current LPM parameter struct (already built from
%                     clinical data via build_patient_params)
%   optimized_table - MATLAB table with columns 'Parameter' and 'Optimized'
%                     as produced by run_lbfgsb_optimization_pda_coa.m
%                     (columns: Parameter, Baseline, Optimized, LB, UB, ChangePercent)
%   disease_mode    - 'PDA_only' or 'PDA_CoA'
%
% OUTPUTS:
%   params       - updated parameter struct with optimized values applied
%   applied_list - cell array of parameter names that were applied
%   ignored_list - cell array of parameter names that were skipped
%
% AUTHOR:   Cardiovascular Simulation Team
% DATE:     2026-06-24
% VERSION:  1.0
% -----------------------------------------------------------------------

applied_list = {};
ignored_list = {};

% Parameters that are only meaningful with an active CoA module
coa_only_params = {'stenosis_pct', 'coa_length_mm'};

% Parameters applied conditionally (PDA_CoA only; useful but not CoA-geometry)
pda_coa_preferred = {'R_pa', 'R_shunt_pda'};

for k = 1:height(optimized_table)

    pname = char(optimized_table.Parameter{k});
    val   = optimized_table.Optimized(k);

    if ~isnumeric(val) || isnan(val)
        continue;
    end

    switch pname

        case 'Emax_lv'
            params.Emax_lv = max(0.01, val);
            params.Emin_lv = params.Emax_lv * 0.05;
            params.Emax_rv = params.Emax_lv * 0.5;
            params.Emin_rv = params.Emax_rv * 0.05;
            applied_list{end+1} = 'Emax_lv'; %#ok<AGROW>

        case 'R_systemic'
            params.R_systemic = max(0.01, val);
            applied_list{end+1} = 'R_systemic'; %#ok<AGROW>

        case 'C_sys'
            params.C_sys = max(1e-4, val);
            applied_list{end+1} = 'C_sys'; %#ok<AGROW>

        case 'C_ao'
            params.C_ao = max(1e-4, val);
            applied_list{end+1} = 'C_ao'; %#ok<AGROW>

        case 'R_pa'
            if strcmp(disease_mode, 'PDA_CoA')
                params.R_pa = max(0.001, val);
                applied_list{end+1} = 'R_pa'; %#ok<AGROW>
            else
                ignored_list{end+1} = pname; %#ok<AGROW>
            end

        case 'R_shunt_pda'
            if strcmp(disease_mode, 'PDA_CoA')
                params.R_shunt_pda = max(0.01, val);
                applied_list{end+1} = 'R_shunt_pda'; %#ok<AGROW>
            else
                ignored_list{end+1} = pname; %#ok<AGROW>
            end

        case 'stenosis_pct'
            if strcmp(disease_mode, 'PDA_CoA')
                params.stenosis_pct = max(0.1, min(99.9, val));
                applied_list{end+1} = 'stenosis_pct'; %#ok<AGROW>
            else
                ignored_list{end+1} = pname; %#ok<AGROW>
            end

        case 'coa_length_mm'
            if strcmp(disease_mode, 'PDA_CoA')
                applied_list{end+1} = 'coa_length_mm'; %#ok<AGROW>
            else
                ignored_list{end+1} = pname; %#ok<AGROW>
            end

        otherwise
            % Unknown parameter — skip without error
    end
end

end
