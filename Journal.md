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