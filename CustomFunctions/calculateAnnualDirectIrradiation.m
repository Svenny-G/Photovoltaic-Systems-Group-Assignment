function G_direct = calculateAnnualDirectIrradiation(s_ix, orientation, G_Bn, Az, hs)
% Computes hourly direct POA irradiance [W/m²] for each module on a roof segment
%
% Inputs:
%   s_ix        - Roof segment index (1–8)
%   orientation - 'portrait' or 'landscape'
%   G_Bn        - 8760x1 beam normal irradiance [W/m²]
%   Az          - 8760x1 solar azimuth [°]
%   hs          - 8760x1 solar altitude [°]
%
% Output:
%   G_direct_hourly - [8760 × N] matrix of direct irradiance per module

    % Determine facing direction and tilt based on segment index
    dirs = ["south", "north", "west", "east"];
    facing_map = [1 2 3 4 1 2 3 4];
    tilt_map = [14.14 14.14 59.29 59.29 14.14 14.14 59.29 59.29];

    facing_direction = dirs(facing_map(s_ix));
    tilt_deg = tilt_map(s_ix);
    az_module_deg = struct('south', 0, 'north', 180, 'east', -90, 'west', 90);
    az_module = az_module_deg.(char(facing_direction));

    % Calculate cos(AOI) for every hour (shared across modules)
    cosAOI = calculateCosAOI(Az, hs, tilt_deg, az_module);
    cosAOI(cosAOI < 0) = 0;

    % Load skylines for this segment and orientation
    file_path = fullfile('MatlabFunctions', 'Building', ...
        sprintf('%s_skylines_%s.mat', orientation, facing_direction));
    data = load(file_path, 'skylines');
    skylines_segment = data.skylines{s_ix};
    n_modules = length(skylines_segment);

    G_direct_module = zeros(8760, n_modules);

    for i = 1:n_modules
        sf = calculateShadingFactor(skylines_segment{i}, Az(:), hs(:))';
        G_direct_module(:, i) = sf .* G_Bn .* cosAOI;
    end
    G_direct = sum(G_direct_module,2);
end
