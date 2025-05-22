function BarPlotMonthavgModuleWorkingTemp(TM,G)
% Inputs: TM (8760x81), G (8760x81)
%
% Makes barplot of average daytime (inferred from nonzero irradience) module temperature for each month
% Credits to ChatGPT for writing most of this function for me
%
% Inputs:
%   TM          - 8760x81 array of calculated module temperature for each module for every hour of the year
%   G           - 8760x81 array of calculated irradience for each module for every hour of the year
%
% Output:
%   barplot of average daytime module temperature per month


% Generate datetime array for the year 2005
t = datetime(2005,1,1,0,0,0) + hours(0:8759);
months = month(t);  % Get month number (1-12) for each hour

avg_monthly_temp = zeros(12,1);

for m = 1:12
    % Indices for current month
    idx = months == m;
    
    % Extract relevant data for this month
    TM_month = TM(idx, :);
    G_month = G(idx, :);
    
    % Mask where irradiance is zero
    mask = G_month > 0;
    
    % Set TM to NaN where irradiance is zero
    TM_month(~mask) = NaN;
    
    % Compute mean temperature excluding NaNs (across time and modules)
    avg_monthly_temp(m) = mean(TM_month(:), 'omitnan');
end

% Plotting
bar(avg_monthly_temp);
xlabel('Month');
ylabel('Average Module Temperature (°C)');
title('Average Monthly Daytime Module Working Temperature');
xticks(1:12);
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
grid on;
end