function ax = plot3DBuildings(vertices,faces,varargin)
% Sintax:
% ax = plot3DBuildings(vertices,faces);
% ax = plot3DBuildings(vertices,faces,'zmin',z1,'zmax',z2,'xmin',x1,'xmax',x2,'ymin',y1,'ymax',y2);

% Description: Plots the 3D model of the building. The axes' limits are set
% to show only the building of interest by default. The axes' limits can be
% modified to show the surrounded buildings.
%
% Examples of use:
% To plot the entire set of building
% load('building2020.mat','building_faces','building_vertices');
% z1 = 0 ; z2 = 50; x1 = -50; x2 = 50; y1 = -50; y2 = 50; 
% ax = plot3DBuildings(building_vertices,building_faces,'zmin',z1,'zmax',z2,'xmin',x1,'xmax',x2,'ymin',y1,'ymax',y2);

% To plot only the building of interest:
% load('building2020.mat','building_faces','building_vertices');
% ax = plot3DBuildings(building_vertices,building_faces);

zmin = 0; zmax = 40; xmin = -9.1; xmax = 7.1; ymin = -3.1; ymax = 35.1;
if ~isempty(varargin)
    try
        for k = 1:2:length(varargin)
            switch(varargin{k})
                case 'xmin'
                    xmin = varargin{k+1};
                case 'xmax'
                    xmax = varargin{k+1};
                case 'ymin'
                    ymin = varargin{k+1};
                case 'ymax'
                    ymax = varargin{k+1};
                case 'zmin'
                    zmin = varargin{k+1};
                case 'zmax'
                    zmax = varargin{k+1};
                otherwise
                    error('Invalid input argument')
            end
            
        end
    catch
        error('Invalid axes limits')
    end
end

fig = figure('units','normalized','outerposition',[0 0 1 1]);
ax = axes(fig);
ax.NextPlot = 'add';

patch('Faces',faces,'Vertices',vertices,'FaceVertexCData', ...
    [0.3 0.3 0.3],'FaceColor','flat','EdgeAlpha',0.15);
axis equal
ax.XLim = [xmin xmax];
ax.YLim = [ymin ymax];
ax.ZLim = [zmin zmax];
end
