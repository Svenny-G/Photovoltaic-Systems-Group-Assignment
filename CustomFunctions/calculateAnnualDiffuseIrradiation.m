function G_diffuse_ann = calculateAnnualDiffuseIrradiation(s_ix, orientation, facing_direction, G_Dh)
% Calculate annual diffuse irradiation per module for a given roof segment and configuration
%
% Inputs:
%   s_ix             - Roof segment index (1 to 8)
%   orientation      - 'portrait' or 'landscape'
%   facing_direction - 'south', 'east', 'west', or 'north'
%   G_Dh             - 8760x1 vector of diffuse horizontal irradiance [W/m²]
%
% Output:
%   G_diffuse_ann    - Nx1 vector of annual diffuse irradiation [kWh/m²] for each module

    % Load sky view factors
    filename = sprintf('%s_skylines_%s.mat', orientation, facing_direction);
    data = load(filename, 'svf');
    svf_segment = data.svf{s_ix};  % cell array, one entry per module

    n_modules = length(svf_segment);
    G_Dh_kWh = G_Dh * 1e-3;  % convert from W to kWh

    G_diffuse_ann = zeros(n_modules, 1);
    for i = 1:n_modules
        % SVF is a constant scalar per module
        svf_i = svf_segment{i};
        % Element-wise multiply and sum over the year
        G_diffuse_ann(i) = sum(svf_i .* G_Dh_kWh);
    end
end
