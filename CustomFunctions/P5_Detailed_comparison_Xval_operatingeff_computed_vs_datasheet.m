%% PROBLEM 5
% Report the average operating module efficiency (including a rough approximation of system
% losses) and the installed PV power (in kilowatts peak) required to generate X% of the electricity
% consumed annually. Justify your approximation of the system losses. This step will give you a general
% idea of the required number of modules for your case.

%plan: first calculate average operating module efficiency for each hour of the year using TC Pmax and efficiency from the datasheet, then calculate expected yield per panel, assume external DC/AC efficiency 
%assume FF = 74% (rough average between STC=75.0% and NOCT=73.3%), assume ideal diode


function P5_Detailed_comparison_Xval_operatingeff_computed_vs_datasheet(monthly_demand,FF,n_modules,Module_temps,G_module_raw,losses_approx,X_val)
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
Mod_eff = Mod_Pmpp./G_mod_mask/Am;

valid_hours = ~isnan(Mod_eff) & G_module_raw > 0 & Mod_eff > 0; 
avg_operating_efficiency = mean(Mod_eff(valid_hours));
avg_annual_irradiation = mean(sum(G_module_raw))*1e-3;
avg_annual_energy_dc = mean(avg_annual_irradiation * avg_operating_efficiency * Am);


%Approximated system losses
losses_approx = losses_approx;        

avg_effective_efficiency = avg_operating_efficiency*(1-losses_approx);
annual_yield = avg_annual_irradiation*avg_effective_efficiency;


%Grid independence levels -> X%
X_val = [10, 20, 25, 30, 45, 60, 65, 80, 100];

LenX = length(X_val);
Req_P = zeros(LenX, 1);
Num_Mod = zeros(LenX, 1);
Req_generation_X = zeros(LenX, 1);
Act_P = zeros(LenX, 1);

for i = 1:LenX
    X = X_val(i);
    Req_generation_X(i) = annual_demand * (X/100);
    Req_P(i) = Req_generation_X(i) / annual_yield;
    Num_Mod(i) = ceil(Req_P(i)*1000/STC_Pmod);
    Act_P(i) = Num_Mod(i) * STC_Pmod / 1000;
end

Q5_resultsA_table = table(X_val', Req_generation_X, Req_P, Num_Mod, Act_P, ...
    'VariableNames', {'Grid independence (%)', 'Req. generation (kWh/year)', ...
                      'Req. power (kWp)', 'Number of modules (#)', 'Actual installed power (kWp)'});

disp('=== RESULTS GRID INDEPENDENCE LEVELS WITH OPERATING EFFICIENCY COMPUTED ===');
disp(Q5_resultsA_table);


%DATASHEET MODULE OPERATING EFFICIENCY------------------------------------- 
valid_hours2 = G_module_raw > 0;
avg_Tmod= mean(T_mod(valid_hours2));
T_derating = 1 + TC_P * (avg_Tmod - STC_T);


%Approximated system losses
losses_approx2 = losses_approx;

avg_operating_efficiency2 = STC_efficiency_mod * T_derating * (1-losses_approx2);
annual_yield2 = avg_annual_irradiation * avg_operating_efficiency2;


%Grid independence levels -> X%
X_val2 = [10, 20, 25, 30, 45, 60, 65, 80, 100];

LenX2 = length(X_val2);
Req_P2 = zeros(LenX2, 1);
Num_Mod2 = zeros(LenX2, 1);
Req_generation_X2 = zeros(LenX2, 1);
Act_P2 = zeros(LenX2, 1);

for i = 1:LenX2
    X = X_val2(i);
    Req_generation_X2(i) = (X / 100) * annual_demand;
    Req_P2(i) = Req_generation_X2(i) / annual_yield2;
    Num_Mod2(i) = ceil(Req_P2(i) * 1000 / STC_Pmod);
    Act_P2(i) = Num_Mod2(i) * STC_Pmod / 1000;
end

Q5_resultsB_table = table(X_val2', Req_generation_X2, Req_P2, Num_Mod2, Act_P2, ...
    'VariableNames', {'Grid independence (%)', 'Req. generation (kWh/year)', ...
                      'Req. power (kWp)', 'Number of modules (#)', 'Actual installed power (kWp)'});

disp('=== RESULTS GRID INDEPENDENCE LEVELS WITH OPERATING EFFICIENCY FROM DATASHEET ===');
disp(Q5_resultsB_table);
end