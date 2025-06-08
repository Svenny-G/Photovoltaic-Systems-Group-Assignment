function [G_total, G_module_matrix, G_module_raw] = calculateTotalIrradiation(s_ix, orientation, G_Bn, G_Dh, G_Gh, Az, hs)
% Computes hourly total POA irradiance [W/m²] for all modules on a roof segment
%
% Inputs:
%   s_ix        - Roof segment index (1–8)
%   orientation - 'portrait' or 'landscape'
%   G_Bn        - 8760x1 beam normal irradiance [W/m²]
%   G_Dh        - 8760x1 diffuse horizontal irradiance [W/m²]
%   G_Gh        - 8760x1 global horizontal irradiance [W/m²]
%   Az          - 8760x1 solar azimuth [°]
%   hs          - 8760x1 solar altitude [°]
%
% Output:
%   G_total     - 8760x1 vector of total irradiance summed over all modules
    Am = 1.7;    % module area
    Az_fix = Az+180;
    % Map segment to facing direction and tilt
    dirs = ["south", "north", "west", "east"];
    facing_map = [1 2 3 4 1 2 3 4];
    tilt_map = [14.14 14.14 59.29 59.29 14.14 14.14 59.29 59.29];

    facing_direction = dirs(facing_map(s_ix));
    tilt = tilt_map(s_ix);
    az_module_deg = struct('south', 0, 'north', 180, 'east', -90, 'west', 90);
    az_module = az_module_deg.(char(facing_direction));

    % Load skylines and SVF
    file_path = fullfile('MatlabFunctions', 'Building', ...
        sprintf('%s_skylines_%s.mat', orientation, facing_direction));
    data = load(file_path, 'skylines', 'svf');
    skylines_segment = data.skylines{s_ix};
    svf_segment = data.svf{s_ix};
    n_modules = length(skylines_segment);

    % Calculate cos(AOI) for all hours (shared)
    cosAOI = calculateCosAOI(Az, hs, tilt, az_module);
    cosAOI(cosAOI < 0) = 0;

    % Allocate matrices
    G_direct_module = zeros(8760, n_modules);
    G_diffuse_module = zeros(8760, n_modules);
    G_albedo_module = zeros(8760, n_modules);

    % Fixed albedo value
    albedo = 0.15;

    for i = 1:n_modules
        % Shading factor for each module
        sf = calculateShadingFactor(skylines_segment{i}, Az(:), hs(:))';

        % Direct
        % direct irradiance = shading factor * beam normal irradiance * cos(incidence angle)
        G_direct_module(:, i) = sf .* G_Bn .* cosAOI;

        % Diffuse
        % diffuse irradiance = sky view factor * diffuse horizontal irradiance
        svf_i = svf_segment(i);  % scalar
        G_diffuse_module(:, i) = svf_i .* G_Dh;

        % Albedo
        % G_albedo = GHI * albedo * (1-SVF)
        G_albedo_module(:, i) = (1 - svf_i) .* albedo .* G_Gh;
    end
    
    % Total annual irradiance per module in n_modules
    G_module_matrix = sum(G_direct_module + G_diffuse_module + G_albedo_module, 1) / Am;  % [1 × N_modules]

    % Total POA irradiance per hour (summed over modules)
    G_total = sum(G_direct_module + G_diffuse_module + G_albedo_module, 2) / Am;

    % Later problems require also hourly data per module
    G_module_raw = (G_direct_module + G_diffuse_module + G_albedo_module) / Am;
end
