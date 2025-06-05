# Journal PV Systems
Use this yournal to document your progress and plans from day to day. From writing/editing functions to progress.

## Problem 1
used provided function GetLoadProfileLatitude.m and location file aberdeen.met to get required data for plotting

## Problem 2
created data files landscape_skylines_'%s'.mat ['north', 'south', 'west'] and portrait_skylines_'%s'.mat ['north', 'south', 'west']
by using the provided function orientation_change.m on the provided data files landscape_skylines_east.mat and portrait_skylines_east.mat

implemented custom functions calculateTotalIrradiation.m and calculateCosAOI.m to calculate the total Irradiance for a chosen roof segment
and solar module orientation. this uses the provided function calculateShadingFactor.m, too.

used provided function plotModulesOnRoof.m to plot module irradiation for all 8 roof segments. Note: couldn't find provided color bar limits
Note2: plot of segment 8 seems weird

## Problem 3
Using the results from problem 2, the Irradience for position 1 in the relevant sector is calculated, from which the expected yearly energy yield is calculated for each of the four models of solar panels. From these results the yearly energy yield per € was calculated for each of the models after which the most cost-effective model was found to be the "Solar Tech TS60-6M3-280S"

## Problem 4
Calculated module temperatures for each module and each hour of the day using the Faiman model and validated the results against the limited data available from the datasheet (temperature is below NOCT for most of the year, which is expected from the weathercoditions in Aberdeen). Used this data to generate a barplot of average working temperature of the modules for each month of the year.

## Problem 5