import Matlabfunctions.*
import CustomFuncs.*

%% PROBLEM 1
%  monthly electricity demand on the left y-axis and the monthly horizontal radiation on the right y-axis.
%  12 months on the x-axis.

% load location data
load('Locations/Aberdeen.mat');
T_Aberdeen = table(N, FF, G_Bn, G_Dh, G_Gh, Az, hs, Ta, Ts);
% disp(T_Aberdeen);
% Variables from Aberdeen.mat
% N:    wtf is this?                        (0-8)
% FF:   Far-field shading factor            (dimensionless)
% G_Bn: Beam Normal Irradiance              (W/m^2)
% G_Dh: Diffuse Horizontal Irradiance       (W/m^2)
% G_Gh: Global Horizontal Irradiance        (W/m^2)
% Az:   solar azimut angle (-180 to 180)    (degrees)                  
% hs:   solar altitude angle/height of sun  (degrees)
% Ta:   ambient temperature                 (° C)
% Ts:   surface temperature                 (° C)
Aberdeen_Latitude = 57.150278;  % from wikimedia – geohack
demand_kWh = GetLoadProfileLatitude(Aberdeen_Latitude) * 1e-3;    % assuming (W), converting to kWh (h=1)

% Generate datetime vector
time = datetime(2024,1,1,0,0,0) + hours(0:8759);
months = month(time)';  % month index 1–12 for each hour

% Monthly demand (kWh)
monthly_demand = accumarray(months, demand_kWh, [12, 1]); % needs column vectors

% Irradiance (convert G_Gh to kWh/m^2)
monthly_GGh = accumarray(months, G_Gh * 1e-3, [12, 1]);

%% PLOTTING
%  monthly electricity demand on the left y-axis and the monthly horizontal radiation on the right y-axis.
%  12 months on the x-axis.
figure;
months = 1:12;

yyaxis left
bar(months, monthly_demand, 0.4);
ylabel('Electricity Demand (kWh)');

yyaxis right
plot(months, monthly_GGh, '-o', 'LineWidth', 2);
ylabel('Global Horizontal Irradiation (kWh/m²)');

xlabel('Month');
title('Monthly Electricity Demand and Irradiation – Aberdeen');
xticks(1:12);
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
legend('Demand', 'GHI', 'Location', 'northwest');
grid on;


%% PROBLEM 2
% Calculate the annual irradiation on each of the allowed mounting places.
% Take into account direct, diffuse, and reflected components of the POA irradiance.
% Choose mounting orientation based on the higher total annual irradiation.

chosen_orientation = strings(1, 8);  % Store best orientation per segment

% loop over all 8 roof segments
for segment = 1:8
    % Compute total irradiation for each orientation
    [G_portrait_total, G_portrait_per_mod] = calculateTotalIrradiation(segment, 'portrait', ...
        G_Bn, G_Dh, G_Gh, Az, hs);
    [G_landscape_total, G_landscape_per_mod] = calculateTotalIrradiation(segment, 'landscape', ...
        G_Bn, G_Dh, G_Gh, Az, hs);

    % Sum total annual irradiation [Wh/m²]
    total_portrait = sum(G_portrait_total);
    total_landscape = sum(G_landscape_total);

    %test
    %fprintf('Segment %d:\n', segment);
    %fprintf('  Portrait:  %.2f kWh/m²\n', total_portrait * 1e-3);
    %fprintf('  Landscape: %.2f kWh/m²\n', total_landscape * 1e-3);

    % Choose orientation
    if total_portrait >= total_landscape
        chosen_orientation(segment) = "portrait";
        G_best = G_portrait_per_mod;
    else
        chosen_orientation(segment) = "landscape";
        G_best = G_landscape_per_mod;
    end

    % Compute annual irradiation per module [kWh/m²]
    annual_kWh_per_m2 = G_best * 1e-3;

    % Plot only the chosen orientation
    modelfile = sprintf('%s_modules.mat', chosen_orientation(segment));
    m_ix = 1:length(annual_kWh_per_m2);

    color_limits = [300 900];  % exceeds min/max of module irradiations
    
    %test
    %fprintf('  Max module: %.2f kWh/m²\n', max(annual_kWh_per_m2));
    %fprintf('  Min module: %.2f kWh/m²\n', min(annual_kWh_per_m2));

    figure;
    plotModulesOnRoof(modelfile, segment, m_ix, 'irradiation', ...
        annual_kWh_per_m2, color_limits);
    title(sprintf('Segment %d - %s (chosen)', segment, chosen_orientation(segment)));
end
