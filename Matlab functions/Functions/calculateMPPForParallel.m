function [pmpp_sys, impp_sys, vmpp_sys, impp_str] = ...
    calculateMPPForParallel(iscs, vocs, FF, impp_stc, isc_stc)
% Description: returns the maximum power point (mpp) parameters (power -
% pmpp_str; current - impp_str; voltage - vmpp_str) and mpp string current
% of series/parallel connected PV modules given their short circuit
% currents (iscs), open circuit voltages (vocs) and module parameters (fill
% factor - FF; maximum power point current at STC conditions - impp_stc;
% short circuit current at STC conditions - isc_stc). Each row in the input
% iscs and vocs represents a series-connected string and are connected in
% parallel to the other rows in the matrices.
%
% Example of use:
% Find the maximum power point parameters for 15 PV modules connected in
% three strings of 5 modules each
% iscs = [4 2 7.4 7.4 7.6; 8 9 10 11 10; 0.1 0.11 0.1 0.1 0.09];
% vocs = [38.5 37.2 39.6 39.6 39.7; 39.8 40 40.2 40.4 40.2; 31.7 31.9 ...
% 31.7 31.7 31.5];
% FF = 0.779; impp_stc = 9.46; isc_stc = 10;
% [pmpp_str, impp_sys, vmpp_str, impp_str] = ...
%   calculateMPPForParallel(iscs, vocs, FF, impp_stc, isc_stc)

% Sort (Voc, Isc) pairs according to decreasing current
[~,idx1] = sort(iscs,2,'descend');
[nn, mm] = size(iscs);
idx1_lin = sub2ind(size(iscs), repmat((1:nn)',1,mm), idx1);
isc_mps = iscs(idx1_lin);
voc_mps = cumsum(vocs(idx1_lin),2);

% Find the pair that gives approximately the maximum power of the system
p_FF_mps = isc_mps.*voc_mps*FF;
[pmpp_sys,idx2] = max(sum(p_FF_mps),[],2);

impp_str = impp_stc*iscs(:,idx2)/isc_stc;
impp_sys = sum(impp_str);
vmpp_sys = pmpp_sys./impp_sys;

end