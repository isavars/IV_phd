function []=lever_remapping_oct2011(varargin)
% Field similarity/remapping, written for Colin Lever data, October 2011.
% Designed to analyse datasets with trial structure:
%
%   A,A,A,B,B,A
%
% Where A is baseline 60cm square, B can be either a 77cm diameter circle
% or a 82cm diameter open platform.
%
% If called with the argument:
%   lever_remapping_oct2011('drawfig')
% it will draw the transformed rate maps used for the analysis.

SD=gss;

% Generate matrix for comparing trials. Basically everything against everything else,
% but slightly complicated by the fact that some datasets will have only 6 trials, some 8, but we want the 
% comparisons for first 6 to always be in the same columns.
count=1;
for ii=1:5
    for jj=ii+1:6
        trialCompMatrix(1:2,count)=[ii;jj]; %#ok<AGROW>
        count=count+1;
    end
end
trialCompLimitShort=length(trialCompMatrix);
trialCompMatrix=[trialCompMatrix,[7 7 7 7 7 7 7 8 8 8 8 8 8 8; 1 2 3 4 5 6 8 1 2 3 4 5 6 7]];
trialCompLimitLong=length(trialCompMatrix);
% Make labels for excel output %
trialLabels = [cellstr(num2str(trialCompMatrix(1,:)')) cellstr(num2str(trialCompMatrix(2,:)'))];
ind=strcmp('7',trialLabels); trialLabels(ind)={'4_scaled'};
ind=strcmp('8',trialLabels); trialLabels(ind)={'5_scaled'};
for ii=1:size(trialLabels,1)
	trialCompLabels{ii} = [trialLabels{ii,1} ' vs ' trialLabels{ii,2}]; %#ok<AGROW>
end
        

% Get user parameters for making rate maps %
prms=rates_params();
if isempty(prms.crop);   prms.crop = [512 512];  prms.abs_space=0;   end    % Force centred crop to align environments.

excelRowCount=2;
for ii=1:length(SD.selData)
    data=evalin('base',SD.selData{ii});
    pearsonCorr=[];
    popVect=[];
    maps=rates_main(data,prms);
    
    % Position Data adjustment for trials in env  'B' %
    templateMapIndex = 3;
    if strcmp(data.trials(4).user.environment,'square')
        % Square data (i.e. no novel probe) - run control topological transform to transform everything to square  %
        for jj=1:6
            maps(jj,:) = rates_transform(maps(jj,:), data.trials(jj).user.shape, maps{templateMapIndex,1}, data.trials(templateMapIndex).user.shape, 'reg_size', 1);
        end
        trialCompLimit=trialCompLimitShort;
        
    elseif strcmp(data.trials(4).user.environment,'circle')
        % Circle data - run topological transform to transform everything to square  %
        for jj=1:6
            maps(jj,:) = rates_transform(maps(jj,:), data.trials(jj).user.shape, maps{templateMapIndex,1}, data.trials(templateMapIndex).user.shape, 'reg_size', 1);
        end
        trialCompLimit=trialCompLimitShort;
        
    elseif strcmp(data.trials(4).user.environment,'open_platform')
        % 82cm platform. Two analyses here: 1) leave data untransformed, 2) shrink by 60%, so create duplicate transformed maps. %
        % Here, we shrink the open platform, and create duplicate trials at the end of the dataset %
        for jj=4:5
            maps(jj+3,:) = rates_transform(maps(jj,:), data.trials(jj).user.shape, maps{templateMapIndex,1}, data.trials(templateMapIndex).user.shape, 'force_scale', 0.6);  % note 'jj+3', to add transformed maps to the end, not replace existing.
        end
        % Here, we 'transform' the open platorm without scaling, i.e. just sample the central portion %
        for jj=4:5
            maps(jj,:) = rates_transform(maps(jj,:), data.trials(jj).user.shape, maps{templateMapIndex,1}, data.trials(templateMapIndex).user.shape, 'force_scale', 1); 
        end
        % Here, we run the 'control transform' on everything else %
        for jj=[1 2 3 6]
            maps(jj,:) = rates_transform(maps(jj,:),data.trials(jj).user.shape,maps{templateMapIndex,1},data.trials(templateMapIndex).user.shape,'reg_size',1);
        end
        trialCompLimit=trialCompLimitLong;
    end
    
    % Draw maps (for testing only ) %
    if ~isempty(varargin) && strcmp(varargin{1},'drawfig')
        data.rate_maps(end+1).maps=maps;
        data.rate_maps(end).name='trans_test';
        gra_mapfig(data, length(data.rate_maps), 1:size(maps,1), 1:size(maps,2));
    end
    
    % Run cell firing comparison analysis % 
    for jj=1:trialCompLimit   
        trA=trialCompMatrix(1,jj);   trB=trialCompMatrix(2,jj);
        
        % Pearsons r correlation of firing rates of bins %
        for kk=1:length(data.trials(1).cells)
            pearsonCorr(kk,jj) = map_spatialcorr( maps{trA,kk}, maps{trB,kk} );
        end    

        % Ensemble population vector comparison %
        popVect(1,jj) = map_popvect( maps(trA,:), maps(trB,:) );   

    end
    
    % Write results to excel %
    % Pearson Corr %
    xlswrite('lever_remapping_oct2011.xlsx', pearsonCorr, 'Pearson r', ['C' num2str(excelRowCount)]);                         % The data
    xlswrite('lever_remapping_oct2011.xlsx', getdatanames(data,'cell','long')', 'Pearson r', ['B' num2str(excelRowCount)]);    % Cell labels
    xlswrite('lever_remapping_oct2011.xlsx', repmat(SD.selData(ii),length(data.trials(1).cells),1), 'Pearson r', ['A' num2str(excelRowCount)]);  % Dataset labels
    excelRowCount=excelRowCount+size(pearsonCorr,1)+1;
    
    % Pop Vect %
    xlswrite('lever_remapping_oct2011.xlsx', popVect, 'Pop Vect', ['C' num2str(ii+1)]);   % The data.
    xlswrite('lever_remapping_oct2011.xlsx', SD.selData(ii), 'Pop Vect', ['A' num2str(ii+1)]);    % The name of the dataset.
    
end

% Label trial comparisons in dataset %
xlswrite('lever_remapping_oct2011.xlsx', trialCompLabels, 'Pearson r', 'C1');
xlswrite('lever_remapping_oct2011.xlsx', trialCompLabels, 'Pop Vect', 'C1');   
















