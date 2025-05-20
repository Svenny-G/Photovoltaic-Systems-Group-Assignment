function cosAOI = calculateCosAOI(Az, hs, tilt, az_module)
% Computes cos(AOI) between sun and module normal for every hour
%
% Inputs:
%   Az            - Solar azimuth [°], 8760x1
%   hs            - Solar elevation [°], 8760x1
%   tilt          - Tilt angle of the module [°]
%   az_module     - Azimuth of module facing direction [°]
%
% Output:
%   cosAOI        - Cosine of angle of incidence, 8760x1

    % Sun vector for each hour
    s_x = cosd(hs) .* sind(Az);  % East-West
    s_y = cosd(hs) .* cosd(Az);  % North-South
    s_z = sind(hs);              % Up

    S = [s_x, s_y, s_z];  % [8760 x 3]

    % Module normal vector (same for all hours)
    n = [sind(tilt) * sind(az_module), ...
         sind(tilt) * cosd(az_module), ...
         cosd(tilt)];  % 1x3

    cosAOI = S * n';
end
