function [pmpp_str, impp_str, vmpp_str] = calculateMPPForSeries(iscs, vocs, ...
    FF, impp_stc, isc_stc)
% Description: returns the maximum power point parameters (power -
% pmpp_str; current - impp_str; voltage - vmpp_str) of series-connected PV
% modules given their short circuit currents (iscs), open circuit voltages
% (vocs) and module parameters (fill factor - FF; maximum power point
% current at STC conditions - impp_stc; short circuit current at STC
% conditions - isc_stc). Each row in the input iscs and vocs represents
% a PV module in the series-connected string. iscs and vocs can have more
% than one row, being each row a different string or time step.
%
% Example of use:
% Find the maximum power point parameters for the example in the lecture
% "Bypass diodes and Mismatch losses"
% iscs = [4	2 7.4 7.6];
% vocs = [38.5, 37.2, 39.6, 39.7];
% FF = 0.779; impp_stc = 9.46; isc_stc = 10;
% [pmpp_str, impp_str, vmpp_str] = calculateMPPForSeries(iscs, vocs, FF, ...
% impp_stc, isc_stc)
%
% Find the maximum power point parameters for 5 series-connected PV modules
% in three different moments in time
% iscs3 = [4 2 7.4 7.4 7.6; 8 9 10 11 10; 0.1 0.11 0.1 0.1 0.09];
% vocs3 = [38.5 37.2 39.6 39.6 39.7; 39.8 40 40.2 40.4 40.2; 31.7 31.9 ...
% 31.7 31.7 31.5];
% FF = 0.779; impp_stc = 9.46; isc_stc = 10;
% [pmpp_str, impp_str, vmpp_str] = calculateMPPForSeries(iscs3, vocs3, FF, ...
% impp_stc, isc_stc)

% Sort (Voc, Isc) pairs according to decreasing current
[~,idx1] = sort(iscs,2,'descend');
[nn, mm] = size(iscs);
idx1_lin = sub2ind(size(iscs), repmat((1:nn)',1,mm), idx1);
isc_mps = iscs(idx1_lin);
voc_mps = cumsum(vocs(idx1_lin),2);

% Find the pair that gives approximately the maximum power of the string
p_noFF_mps = isc_mps.*voc_mps;
[~,idx2] = max(p_noFF_mps,[],2);
idx2_lin = sub2ind(size(iscs), (1:nn)', idx2);
idx3 = idx1(idx2_lin);
idx3_lin = sub2ind(size(iscs), (1:nn)', idx3);

% Find the pmpp parameters
pmpp_str = FF*p_noFF_mps(idx2_lin);
impp_str = impp_stc*iscs(idx3_lin)/isc_stc;
vmpp_str = pmpp_str./impp_str;

end