function plotModulesOnRoof(filename,s_ix,m_ix,varargin)
% Sintax:
% plotModulesOnRoof(filename,s_ix,m_ix)
% plotModulesOnRoof(filename,s_ix,m_ix,'irradiation',x)
% plotModulesOnRoof(filename,s_ix,m_ix,'irradiation',x,ls)
% plotModulesOnRoof(filename,s_ix,m_ix,'svf',x)
% plotModulesOnRoof(filename,s_ix,m_ix,'svf',x,ls)
%
% Description:
% This function plots the PV modules on one roof segment.
% filename can be either 'portrait_modules.mat' or 'landscape_modules.mat'
% s_ix is an index (from 1 to 8) that indicates the roof segment and m_ix
% is a vector that indicates the index of the desired modules (see PPT).
% By default the function plots the modules in blue.
% x is a vector with the annual irradiation on every module expressed in
% MWh/m2 or a vector with the sky view fator of every module. The vector
% must be of the same length as m_ix.
% ls can be used to set the the maximum and minimum limits of the colorbar.
% To obtain consistent colors, it is recommended to use ls when this
% function is called multiple consecutive times, e.g. when plotting modules 
% in different roof sectors. 
%
% Example of use 1: Load and plot the building model and then draw modules
% 4 7 and 12 on roof segment 5 in landscape orientation. The annual
% irradiation on the modules is 0.7 0.8 and 1 MWh/m2, respectively.
%
% load('building2020.mat','building_faces','building_vertices');
% ax = plot3DBuildings(building_vertices,building_faces);
%
% s_ix = 5;
% m_ix = [4 7 12];
% irrs = [0.7 0.8 1];
% plotModulesOnRoof('landscape_modules',s_ix,m_ix,'irradiation',irrs);
%
% Example of use 2: Load and plot the building model and then draw modules
% 1 and 2 on roof segment 1 in landscape orientation and module 10 on roof
% segment 2 in portrait orientation.
% irradiation on the modules is 0.5 0.8 and 1.1 MWh/m2, respectively.
% The limits for the colorbar are 0.2 and 1.5 MWh/m2.
%
% load('building2020.mat','building_faces','building_vertices');
% plot3DBuildings(building_vertices,building_faces);
% cb_limits = [0.2 1.5];
% s_ix = 1;
% m_ix = [1 2];
% irrs = [0.5 0.8];
% plotModulesOnRoof('landscape_modules',s_ix,m_ix,'irradiation',irrs,cb_limits);
% s_ix = 2;
% m_ix = [10];
% irrs = [1.1];
% plotModulesOnRoof('portrait_modules',s_ix,m_ix,'irradiation',irrs,cb_limits);


try
    load(filename,'vnorm','vpoints');
catch
    error(['The file ',filename,' is not in the path']);
end
plotting_vals = false;
nmod = length(m_ix);
if isempty(varargin)
    c_ix = ones(nmod,1);
    colors = [0 0 1];%blue
elseif length(varargin)==2 && (strcmp(varargin{1},'svf') || ...
        strcmp(varargin{1},'irradiation')) && length(varargin{2})==nmod
    values = varargin{2};
    if size(values,1)==1
        values=values'; %make it always a column vector
    end
    NCOLORS = 128;
    colors = parula(NCOLORS);
    if length(values)>1
        min_val = min(values);
    else
        min_val = 0;
    end
    max_val = max(values);
    [~,c_ix] = min(abs(linspace(min_val,max_val,NCOLORS)-values),[],2);
    plotting_vals = true;
elseif length(varargin)==3 && (strcmp(varargin{1},'svf') || ...
        strcmp(varargin{1},'irradiation')) && length(varargin{2})==nmod ...
        && numel(varargin{3})==2
    values = varargin{2};
    if size(values,1)==1
        values=values'; %make it always a column vector
    end
    ls = varargin{3};
    NCOLORS = 128;
    colors = parula(NCOLORS);
    min_val = ls(1);
    max_val = ls(2);
    [~,c_ix] = min(abs(linspace(min_val,max_val,NCOLORS)-values),[],2);
    plotting_vals = true;
else
    error('Aditional parameters are invalid');
end

switch filename
    case {'portrait_modules','portrait_modules.mat'}
        R = [-0.5 1.7/2 0
            0.5 1.7/2 0
            0.5 -1.7/2 0
            -0.5 -1.7/2 0]';
    case {'landscape_modules','landscape_modules.mat'}
        R = [-1.7/2 0.5 0
            1.7/2 0.5 0
            1.7/2 -0.5 0
            -1.7/2 -0.5 0]';
    otherwise
        error('Invalid filename');
end
ax = gca;
ax.DataAspectRatio = [1 1 1];
ax.FontSize = 16;
ax.Box  ='on';
ax.XTick = [];
ax.XLabel.String = '';
ax.YLabel.String = '';
ax.YTick = [];
ax.ZTick = [];
ax.FontName = 'Calibri';
view(ax,0,90);

for ix=1:nmod
    m = m_ix(ix);
    center = vpoints{s_ix}(m,:)';
    center(3) = center(3)+0.2;
    tilt = acosd(vnorm{s_ix}(m,3));
    azim = atan2d(vnorm{s_ix}(m,1),vnorm{s_ix}(m,2));
    azim = azim+(azim<0)*360;
    R2 = rotz(-azim)*(rotx(-tilt)*R)+center;
    patch(ax,'Vertices',R2','Faces',[1 2 3 4],'FaceColor',colors(c_ix(ix),:),'FaceAlpha',1,'CDataMapping','direct');
end

if plotting_vals
    NTICK = 3;
    cb = colorbar(ax,'Location','southoutside');
    ticks = linspace(cb.Limits(1),cb.Limits(2),NTICK);
    ticks_l = linspace(min_val,max_val,NTICK);
    cb.Ticks = ticks;
    cb.TickLabels = {sprintf('%.2f',ticks_l(1)),sprintf('%.2f',ticks_l(2)),...
        sprintf('%.2f',ticks_l(3))};
    if strcmp(varargin{1},'irradiation')
        cb.Label.String = 'Annual Irradiation (MWh m^{-2} yr^{-1})';
    elseif  strcmp(varargin{1},'svf')
        cb.Label.String = 'Sky View Factor (-)';
    else
        error('Invalid string!');
    end
    cb.Label.FontName = 'Calibri';
    cb.Label.FontWeight = 'Bold';
end

end

%%
function r = rotx(angle)
r = [1 0 0
    0 cosd(angle) -sind(angle)
    0 sind(angle) cosd(angle)];
end
%%
function r = rotz(angle)
r = [cosd(angle) -sind(angle) 0
    sind(angle) cosd(angle) 0
    0 0 1];
end


