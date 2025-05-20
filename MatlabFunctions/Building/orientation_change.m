desired_orientation = 'west';

if strcmp(desired_orientation,'east')
    phi = 0;
elseif strcmp(desired_orientation,'south')
    phi = 270;
elseif strcmp(desired_orientation,'west')
    phi = 180;
elseif strcmp(desired_orientation,'north')
    phi = 90;
else
    error('Unknown desired orientation');
end

for module_orientation = {'landscape','portrait'}
    % Load from east-facing base files
    load(strcat(module_orientation{1},'_skylines_east.mat'),'skylines','svf');

    skylines_2 = skylines;
    for ii = 1:length(skylines)
        skylines_2{ii} = cellfun(@(x) x(:,[(phi+1):end 1:phi]),skylines{ii},...
            'UniformOutput',false);
    end

    skylines = skylines_2;

    % Save rotated version
    save(strcat(module_orientation{1},'_skylines_', desired_orientation, '.mat'),...
        'skylines','svf');
end
