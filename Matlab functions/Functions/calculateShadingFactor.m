function sf = calculateShadingFactor(skyline,sun_azim,sun_alt)
% Sintax:
% sf = calculateShadingFactor(skyline,As,as)
%
% Description: returns the shading factor for the specified skyline profile
% and the specified solar azimuth and altitude values. sun_azim and sun_alt
% can be 1-D vectors so you can calculate the shading factor for an entire 
% year with only one call to this function.
%
% Example of use:
% Calculate if the sun is blocked from the perspective of the first module
% on the first roof sector when the sun As=180 and as=39 and also when 
% As=170 and as=20;
% load('portrait_skylines.mat')
% sf = calculateShadingFactor(skylines{1}{1},[180 170],[39 20])

if numel(sun_alt) ~= numel(sun_azim)
    error('Solar azimuth and solar altitude vectors are not the same size')
end

[ROWS,COLS] = size(skyline);
center_azim = single(linspace(0+360/(2*COLS),360-360/(2*COLS),COLS));
center_alt = single(linspace(90-180/(2*ROWS),-90+180/(2*ROWS),ROWS))';
sf=false(1,numel(sun_alt));
for k=1:numel(sf)
    [~,row] = min(abs(center_alt-sun_alt(k)));
    [~,col] = min(abs(center_azim-sun_azim(k)));
    sf(k) = skyline(row,col);
end

end
