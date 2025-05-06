desired_orientation = 'east';

if strcmp(desired_orientation,'east')
    phi = 90;
elseif strcmp(desired_orientation,'south')
    phi = 180;
elseif strcmp(desired_orientation,'west')
    phi = 270;
end

skylines_2 = skylines;
for module_orientation = {'landscape','portrait'}
    load(strcat(module_orientation{1},'_skylines.mat'),'skylines','svf');
    for ii = 1:length(skylines)
        skylines_2{ii} = cellfun(@(x) x(:,[(phi+1):end 1:phi]),skylines{ii},...
            'UniformOutput',false);
    end
    skylines = skylines_2;
    save(strcat(module_orientation{1},'_skylines_',desired_orientation,'.mat'),...
        'skylines','svf');
end
