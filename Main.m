import MatlabFunctions.*
import CustomFunctions.*

%% PROBLEM 1
%  monthly electricity demand on the left y-axis and the monthly horizontal radiation on the right y-axis.
%  12 months on the x-axis.

% load location data
load('Locations/Aberdeen.mat');
T_Aberdeen = table(N, FF, G_Bn, G_Dh, G_Gh, Az, hs, Ta, Ts);

% disp(T_Aberdeen);
% Variables from Aberdeen.mat
% N:    Cloud cover                         (0-8)
% FF:   Wind speed                          (m/s)
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

Am = 1.7;  % module area

chosen_orientation = strings(1, 8);  % Store best orientation per segment

% loop over all 8 roof segments -> oops, just segment 4 relevant
for segment = 4
    % Compute total irradiation for each orientation
    [G_portrait_total, G_portrait_per_mod, G_module_raw] = calculateTotalIrradiation(segment, 'portrait', ...
        G_Bn, G_Dh, G_Gh, Az, hs);
    [G_landscape_total, G_landscape_per_mod, G_module_raw] = calculateTotalIrradiation(segment, 'landscape', ...
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

    % Compute annual irradiation per module [kWh/m²/y]
    annual_MWh_per_module = G_best * 1e-6 / Am;    % divide by module area

    % Plot only the chosen orientation
    modelfile = sprintf('%s_modules.mat', chosen_orientation(segment));
    m_ix = 1:length(annual_MWh_per_module);

    color_limits = [0.26 0.48];  % exceeds min/max of module irradiations, given in assignment
    
    %test
    fprintf('  Max module: %.2f kWh/m²/y\n', max(annual_MWh_per_module));
    fprintf('  Min module: %.2f kWh/m²/y\n', min(annual_MWh_per_module));

    figure;
    plotModulesOnRoof(modelfile, segment, m_ix, 'irradiation', ...
        annual_MWh_per_module, color_limits);
    title(sprintf('Segment %d - %s (chosen)', segment, chosen_orientation(segment)));
end

%% PROBLEM 3
%calculate the expected yearly energy output for all four PV modules installed in module position 1 of your roof sector.
[G_landscape_total, G_landscape_per_mod, G_module_raw] = calculateTotalIrradiation(4, 'landscape', ...
        G_Bn, G_Dh, G_Gh, Az, hs);
E = G_landscape_per_mod(1)*1.7e-3*[0.172,0.184,0.196,0.19]; %Calculate energy per annum for each of the four models in kwh/y
'module yield per cost estimate'
E./[0.41*280,0.46*310,0.62*330,0.56*310]

areas = [1.627, 1.683, 1.682, 1.658];  % module areas from datasheets, corresponding to correct effs
effs = [0.172, 0.184, 0.196, 0.190];

E_acc_area = G_landscape_per_mod(1) * 1e-3 .* (areas .* effs);  % now per-model
'module yield per cost estimate accounting for module area'
E_acc_area./[0.41*280,0.46*310,0.62*330,0.56*310]
% doesn't change anything wrt model choice, but is probably more precise this way

%% PROBLEM 4
% Knowing the irradiance incident on each PV module, weather parameters and the module of
% your choice, calculate the module temperature for every hour of the year for each of the allowed
% mounting places. Make one bar plot with the average module temperature during the daytime (when
% the modules are working) and the months on the horizontal axis. Discuss how these temperatures will
% affect the power generated by the solar panels compared to the case when the modules operate at a
% constant temperature of 25°C.
[G_landscape_total, G_landscape_per_mod, G_module_raw] = calculateTotalIrradiation(4, 'landscape', ...
        G_Bn, G_Dh, G_Gh, Az, hs);

    % Load sky view factors
    filename = sprintf('%s_skylines_%s.mat', 'landscape', 'east');
    data = load(filename, 'svf');
    svf_segment = data.svf{4};  % cell array, one entry per module

    n_modules = length(svf_segment);
G_module_raw;
%using Faiman model since it is the best performing simple model for non-arid climates
Module_temps = repmat(Ta,1,n_modules)+G_module_raw./(25*ones(length(FF),n_modules)+6.84*repmat(FF,1,n_modules)); %U0 and U1 were assumed to be 4.3 and 1.55 respectively based on values found here: https://pvpmc.sandia.gov/modeling-guide/2-dc-module-iv/module-temperature/faiman-module-temperature-model/
%These parameters also replicate the NOCT cell temperature within 0.4 degrees
BarPlotMonthavgModuleWorkingTemp(Module_temps,G_module_raw);

%% PROBLEM 5
% Report the average operating module efficiency (including a rough approximation of system
% losses) and the installed PV power (in kilowatts peak) required to generate X% of the electricity
% consumed annually. Justify your approximation of the system losses. This step will give you a general
% idea of the required number of modules for your case.

%plan: first calculate average operating module efficiency for each hour of the year using TC Pmax and efficiency from the datasheet, then calculate expected yield per panel, assume external DC/AC efficiency 
%assume FF = 74% (rough average between STC=75.0% and NOCT=73.3%), assume ideal diode

%Datasheet & constants
Am = 1.7;
STC_Pmod = 280;
STC_Voc = 39.56;
STC_Vmp = 31.90;
STC_Isc = 9.46;
STC_Imp = 8.80;
STC_efficiency_mod = 0.1721;
STC_T =25;
TC_Isc = 0.069/100;
TC_Voc = -0.312/100;
TC_P = -0.432/100;             
kb_T = 298.15*1.381e-23;
q = 1.602e-19;

annual_demand = sum(monthly_demand);

%COMPUTED MODULE OPERATING EFFICIENCY-------------------------------------- 
T25 = 25*ones(length(FF),n_modules);

T_mod = Module_temps;

G_mod_mask = G_module_raw;
G_mod_mask(G_mod_mask<=0)=NaN;

Mod_Voc = STC_Voc*ones(length(FF),n_modules)+kb_T/q*T_mod.*log(G_mod_mask/1000)+STC_Voc*TC_Voc*(T_mod-T25);
Mod_Isc = STC_Isc*G_mod_mask/1000 + STC_Isc*TC_Isc*(T_mod-T25);
Mod_Pmpp= 0.74*Mod_Voc.*Mod_Isc;
Mod_Pmpp(Mod_Pmpp<=0)=NaN;
Mod_eff = Mod_Pmpp./G_mod_mask;

valid_hours = ~isnan(Mod_eff) & G_module_raw > 0 & Mod_eff > 0 ;
avg_operating_efficiency = mean(Mod_eff(valid_hours));
% avg_annual_irradiation = mean(sum(G_module_raw))*1e-3;
% avg_annual_energy_dc = mean(avg_annual_irradiation * avg_operating_efficiency * Am); %invalid calculation method: efficiency and irracience are not independent


%Approximated system losses
losses_approx = 0.12;       %decreased losses to 12% as this is more realistic for grid tied PV , based on inverter efficiency of 92% and other losses 4%

%effective_efficiency = mod_eff*(1-losses_approx); redundant
%annual_yield = avg_annual_irradiation*avg_effective_efficiency; %invalid calculation method, see above
annual_yield = (1-losses_approx)*sum(Mod_Pmpp,1,"omitnan")

%Grid independence levels -> X% = 45%
X_45 = 45;

Req_generation_X = annual_demand * (X_45/100)*1.1;
%Num_Mod = Req_generation_X / mean(annual_yield/1000)
%Num_Mod = ceil(Req_P*1000/STC_Pmod); redundant


[selected_modules,num_panels_req,expected_yield] = Panelselection(annual_yield/1000,Req_generation_X)
Req_P = num_panels_req * STC_Pmod / 1000;


Q5_resultsA_table = table(X_45', Req_generation_X, Req_P, num_panels_req, expected_yield, ...
    'VariableNames', {'Grid independence (%)', 'Req. generation (kWh/year)', ...
                      'Installed power (kWp)', 'Number of modules (#)', 'Expeced yield (kWh/year)'});

disp('=== RESULTS GRID INDEPENDENCE LEVELS WITH OPERATING EFFICIENCY COMPUTED ===');
disp(Q5_resultsA_table);

%% P5 OPTIONAL -> uncomment (Comparison of operating efficiency computed vs datasheet)

% X_val = [10, 20, 25, 30, 45, 60, 65, 80, 100];
% losses_approx = 0.15;
% P5_Detailed_comparison_Xval_operatingeff_computed_vs_datasheet(monthly_demand,FF,n_modules,Module_temps,G_module_raw, losses_approx, X_val);

%% PROBLEM 6
% Optimize the distribution and interconnection of PV modules on the rooftop. Explain clearly in
% the report (with a simple schematic) how many arrays (groups) of modules compose your PV system
% and how are the modules in each array connected with each other and to the inverter. Try to minimize
% the losses due to current mismatch between series-connected modules. Use the symbol below to
% represent a PV module. Provide moreover a plot of the annual irradiation of only the modules that
% compose your PV system, in the chosen mounting orientation. The limits of the color bar range will be
% provided and depend on the city.

% === INPUTS from Problem 5 (assumed already in workspace) ===
% G_module_raw          : [8760 × n_modules] POA irradiance [W/m²]
% Mod_Voc               : [8760 × n_modules] module Voc per hour
% STC_Pmod              : Module nominal power (W)
% Am                    : Module area (m²)
% STC_efficiency_mod    : Efficiency at STC
% annual_demand         : Annual electricity consumption [kWh]
% annual_yield          : Per-module expected yield (Wh)
% losses_approx         : Combined system losses fraction
% Panelselection        : Selection function (provided)
% X_45                  : Grid independence level (%)

% select highest annual P modules
Req_generation_X = annual_demand * (X_45 / 100) * 1.1;    % include 10% safety margin (factor 1.1)
[selected_modules, num_panels_req, expected_yield] = Panelselection(annual_yield/1000, Req_generation_X);

% annual irradiation per module [kWh/m²]
annual_irradiation = sum(G_module_raw, 1) * 1e-3;
selected_irradiation = annual_irradiation(selected_modules);

% voltage limit check and string calculation
Voc_max = max(Mod_Voc, [], 'all', 'omitnan');
Max_mod_per_string = floor(1000 / Voc_max * 0.9);  % 1000 V grid connected inverter voltage maximum; include 10% safety margin (factor 0.9)
N_strings = ceil(num_panels_req / Max_mod_per_string);

% assign modules to strings while minimizing losses
sorted_modules = selected_modules;  % already sorted best-to-worst

string_modules = cell(N_strings, 1);
string_irradiation = cell(N_strings, 1);

for s = 1:N_strings
    start_idx = (s - 1) * Max_mod_per_string + 1;
    end_idx = min(s * Max_mod_per_string, num_panels_req);
    
    string_modules{s} = sorted_modules(start_idx:end_idx);
    string_irradiation{s} = selected_irradiation(start_idx:end_idx);
end

% compute mismatch losses
% Assumption: Energy yield is proportional to irradiation * module rated power (P_STC / 1000)
% Assumes all modules operate ideally at STC efficiency (no temperature or system derating here)

% Total ideal energy if each module operates at its own irradiation level
ideal_energy = sum(selected_irradiation) * STC_Pmod / 1000;  % [kWh]

% Total actual energy if each string limited by its lowest-performing module
strings_energy = 0;
for s = 1:N_strings
    min_I = min(string_irradiation{s});     % worst irradiation in the string [kWh/m²]
    n_mod = length(string_modules{s});      % number of modules in the string
    strings_energy = strings_energy + min_I * n_mod * STC_Pmod / 1000;  % [kWh]
end

% Mismatch loss due to series connection (in % of ideal)
mismatch_loss_pct = (ideal_energy - strings_energy) / ideal_energy * 100;


% Summary Table
Q6_summary_table = table(num_panels_req, num_panels_req * STC_Pmod / 1000, N_strings, ...
     strings_energy, ...
    strings_energy / annual_demand * 100, mismatch_loss_pct, ...
    'VariableNames', {'Total modules', 'System power (kWp)', 'Number of strings', ...
    'Annual generation (kWh)', 'Grid independence (excludes losses) (%)', 'Mismatch loss (%)'});
disp('--- SUMMARY TABLE (Problem 6) ---');
disp(Q6_summary_table);

% Breakdown per string
String_ID = (1:N_strings)';
Mean_Irr = zeros(N_strings, 1);
Std_Irr = zeros(N_strings, 1);
Module_List = strings(N_strings, 1);

for s = 1:N_strings
    Mean_Irr(s) = mean(string_irradiation{s});
    Std_Irr(s) = std(string_irradiation{s});
    Module_List(s) = join(string(string_modules{s}), ', ');
end

Q6_strng_table = table(String_ID, Module_List, Mean_Irr, Std_Irr, ...
    'VariableNames', {'String_ID', 'Modules', 'Mean_Irradiation (kWh/m²)', 'Std Dev'});
disp('--- STRING ASSIGNMENT TABLE ---');
disp('Connect in this order (highest to lowest irradiation) to minimize mismatch and maximize energy harvest');
disp(Q6_strng_table);

% Plot all selected modules
color_limits = [0.26 0.48];
modelfile = 'landscape_modules.mat';     % repeat cause I looped over all 8 roof segments earlier

selected_irradiation_per_mod = selected_irradiation * 1e-3 / Am;   % scaled per module and to MWh
for s = 1:N_strings
    string_irradiation_per_mod{s} = string_irradiation{s} * 1e-3 / Am;
end

figure;
plotModulesOnRoof(modelfile, 4, selected_modules, 'irradiation', ...
    selected_irradiation_per_mod, color_limits);
title(sprintf('Segment 4 – Selected %d Modules (%.1f kWp)', ...
    num_panels_req, num_panels_req * STC_Pmod / 1000));
colorbar;
saveas(gcf, fullfile('Figures', 'Problem6.fig'));

% Plot each string individually (not asked by problem statement)
for s = 1:N_strings
    figure('Name', sprintf('String %d', s), 'NumberTitle', 'off');
    plotModulesOnRoof(modelfile, 4, string_modules{s}, 'irradiation', ...
        string_irradiation_per_mod{s}, color_limits);
    title(sprintf('String %d – %d Modules (%.1f ± %.1f kWh/m²)', ...
        s, length(string_modules{s}), Mean_Irr(s), Std_Irr(s)));
end

%% PROBLEM 7
% Keeping in mind the electrical layout of your PV array(s), calculate the
% DC electrical power generated by the solar panels of your PV array(s) for
% every hour during an entire year using climate data and taking into
% account efficiency losses due to temperature effects. Account for the
% mismatch between PV modules by using the code provided6. We want you to
% report the IV curve of each string of your system at 12:00 on 21st June
% and to complete the following table (one row per string).

%% First, calculate total annual DC power generated from selected modules
% under consideration of temperature effects limiting efficiency and
% current mismatch

% Initialize output matrices
Pmpp_hourly = zeros(8760, N_strings);
Impp_hourly = zeros(8760, N_strings);
Vmpp_hourly = zeros(8760, N_strings);

% Calculate Fill Factor from solar module data sheet values
Fill_Factor = (STC_Vmp * STC_Imp) / (STC_Voc * STC_Isc);  % ≈ 0.796


for s = 1:N_strings
    % Extract module indices for the current string
    mod_ids = string_modules{s};

    % Get hourly Isc and Voc for all modules in this string
    iscs = Mod_Isc(:, mod_ids);   % [8760 × n_mod]
    vocs = Mod_Voc(:, mod_ids);   % [8760 × n_mod]

    % Replace NaNs with zeros (non-generating or invalid modules)
    vocs(isnan(vocs)) = 0;
    iscs(isnan(iscs)) = 0;

    % Compute MPP values using provided function
    [Pmpp_str, Impp_str, Vmpp_str] = calculateMPPForSeries(iscs, vocs, ...
        Fill_Factor, STC_Imp, STC_Isc);

    % Store results for this string
    Pmpp_hourly(:, s) = Pmpp_str;  % [W]
    Impp_hourly(:, s) = Impp_str;
    Vmpp_hourly(:, s) = Vmpp_str;
end

% Total generated power
total_Pmpp = sum(Pmpp_hourly, 2);                   % [W]
total_energy_generated = sum(total_Pmpp) / 1000;    % [kWh]

% summary output
fprintf('=== Annual DC Output Summary ===\n');
fprintf('Total DC energy output: %.1f kWh\n', total_energy_generated);
fprintf('DC Grid independence: %.1f %%\n', 100 * total_energy_generated / annual_demand);

%% Second, generate the IV curve of each string at 12:00 noon on June 21st (summer solstice)

% Time index for June 21st, 12:00
dt = datetime(2023, 1, 1, 0, 0, 0) + hours(0:length(Mod_Isc)-1);
target_time = datetime(2023, 6, 21, 12, 0, 0);
time_idx = find(dt == target_time);

for s = 1:N_strings
    module_idx = string_modules{s};

    vocs_now = Mod_Voc(time_idx, module_idx);
    iscs_now = Mod_Isc(time_idx, module_idx);

    % Sort by descending Isc (bypass logic)
    [iscs_sorted, idx] = sort(iscs_now, 'descend');
    vocs_sorted = vocs_now(idx);

    % Compute cumulative voltage and step currents
    Vk = cumsum(vocs_sorted);
    Ik = iscs_sorted;

    % Table to summarize string stats
    IV_table = table('Size', [N_strings, 7], ...
    'VariableTypes', repmat({'double'}, 1, 7), ...
    'VariableNames', {'String', 'ISC_min', 'ISC_max', 'VOC_min', 'VOC_max', 'Pmin', 'Pmax'});

    for s = 1:N_strings
        module_idx = string_modules{s};
        vocs_now = Mod_Voc(time_idx, module_idx);
        iscs_now = Mod_Isc(time_idx, module_idx);
    
        valid = ~isnan(vocs_now) & ~isnan(iscs_now);
        vocs_now = vocs_now(valid);
        iscs_now = iscs_now(valid);
    
        [iscs_sorted, idx] = sort(iscs_now, 'descend');
        vocs_sorted = vocs_now(idx);
        Ik = iscs_sorted;
        Vk = cumsum(vocs_sorted);
        Pk = Ik .* Vk;
        
        IV_table.String(s) = s;
        IV_table.ISC_min(s) = min(Ik);
        IV_table.ISC_max(s) = max(Ik);
        IV_table.VOC_min(s) = min(vocs_sorted);
        IV_table.VOC_max(s) = sum(vocs_sorted);  % full string voltage
        IV_table.Pmin(s) = Fill_Factor * min(Pk);
        IV_table.Pmax(s) = Fill_Factor * max(Pk);
    end

disp('=== IV Table for 21st June 12:00 ===');
disp(IV_table);


    % Add initial and final steps for full stair-step plot
    V_plot = [];
    I_plot = [];

    V_prev = 0;
    for k = 1:length(Vk)
        V_now = Vk(k);
        I_now = Ik(k);
        V_plot = [V_plot, V_prev, V_now];
        I_plot = [I_plot, I_now, I_now];
        V_prev = V_now;
    end

    % Plot stair-step IV curve
    fig = figure('Name', sprintf('String %d IV Stair-Step', s), ...
                 'NumberTitle', 'off');
    plot(V_plot, I_plot, 'LineWidth', 2);
    xlabel('Voltage [V]');
    ylabel('Current [A]');
    title(sprintf('IV Curve, String %d (21st June, 12:00)', s));
    grid on;
    xlim([0, Vk(end)]);
    ylim([Ik(end) * 0.95, Ik(1) * 1.05]);
    legend('IV Curve', 'Location', 'northeast');

    % Save as .fig
    saveas(fig, fullfile('Figures', sprintf('Problem7_String%d.fig', s)));
end

%% PROBLEM 8