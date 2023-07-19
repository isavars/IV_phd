function [Res] = shuffled_rate_maps_by_cell(varargin)
% Spike shuffled thresholds for rate for single rate maps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prms.shufMinOffset = 30;       % In sec
prms.shufNSteps = 100;        % changes to 100 shuffles and 30seconds form knierim 
prms.dirBin = 6;               %
prms.dirSmooth = 5;            %
prms.timeSliceWindow = 60;    % In sec
prms.maxNTrial = 6;
% maps - place params addrd from getSpatData 
prms.placeBin     = 10; %cam pix/spatial bin
prms.placeSmooth  = 5;  % size smoothing kernel (in bins)
prms.adapt_smooth = 200; % 
prms.crop         = [450 450];
prms.speedFilt    = [2.5 400]; % lower/uppper limits cm/s

if ~isempty(varargin)
    for ii=1:(length(varargin)-1)
        prms.(varargin{ii}) = varargin{ii+1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SD=gss;
% rateMapParams = rates_params('bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth,'space','dir');
% posMapParams = rates_params('bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth,'space','dir','mode','pos');
rateMapParams = rates_params('bin',prms.placeBin,'smooth',prms.placeSmooth,'filt_speed',prms.speedFilt,'adaptive_smooth',prms.adapt_smooth,'crop',prms.crop);
posMapParams = rates_params('bin',prms.placeBin,'smooth',prms.placeSmooth,'filt_speed',prms.speedFilt,'crop',prms.crop,'mode', 'pos');
    
for ii=1:length(SD.selData)
    bigDataArray{ii,1} = evalin('base',SD.selData{ii});
    bigDataArray{ii,2} = SD.selData{ii};
end

[rMapTemp,siTemp,stabTemp,datasetNamesTemp] = deal(cell(1,length(SD.selData)));   % Pre-allocate temp cell arrays for results (cannot combine into structure until after PARFOR).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop for getting shuffled data scores %
tic
parfor ii=1:length(SD.selData)
    data=bigDataArray{ii,1};
    
    % Unfortunately, to take account of different trial lengths, we have to spilt the analyis by trial %   
    [siTempDS,stabTempDS] = deal(nan( length(data.trials), length(data.trials(1).cells), prms.shufNSteps));   % Pre-allocate temp arrays for this dataset.
    rMapTempDS = cell(length(data.trials), length(data.trials(1).cells), prms.shufNSteps);
    for jj=1:length(data.trials)  %% jj = iterator for trial (but note that the first action here is to reduce the data to one trial only, so after that jj only used for output assignment).
        trData=data;
        trData.trials = trData.trials(jj);
        offsets = linspace( prms.shufMinOffset, trData.trials.dur-prms.shufMinOffset, prms.shufNSteps);
        
        % Also need to loop through *Each indivdual offset* %
        for kk=1:length(offsets) %% kk = iterator for offset
            
            dataShuf = spk_randomise(trData,'fixedWrap',offsets(kk));
            
            % Make maps (just for 1 trial, one offset %
            rateMaps = rates_main(dataShuf,rateMapParams)
            posMaps = rates_main(dataShuf,posMapParams);
            tempParams = rateMapParams;   tempParams.filt_1stHalf = 1;
            maps1stHalf = rates_main(dataShuf,tempParams);
            tempParams.filt_1stHalf = 0;   tempParams.filt_2ndHalf = 1;
            maps2ndHalf = rates_main(dataShuf,tempParams);
            
            
            % Get the scores for each cell
            for mm=1:length(trData.trials.cells)
%                 dMap = rateMaps{1,mm};
%                 rvTempDS(jj,mm,kk) = dir_rayleighvector(dMap);
                rMap = rateMaps{1,mm};
                %rMapTempDS{jj,mm,kk} = rMap;
                siTempDS(jj,mm,kk) = map_skaggsinfo(rMap,posMaps{1});
                stabTempDS(jj,mm,kk) = map_spatialcorr( maps1stHalf{1,mm}, maps2ndHalf{1,mm} );
            end

        end
        
        % Transfer the temp arrays for this dataset to the overall cell array %
        %rvTemp{ii} = rvTempDS;
        %rMapTemp{ii} =rMapTempDS;
        siTemp{ii} = siTempDS;
        stabTemp{ii} = stabTempDS;
        datasetNamesTemp{ii} = bigDataArray{ii,2};
        
    end
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-organise the output as a table, for ease of comparison with main data %
% Set up results table %
scoreDum = cell(1,prms.maxNTrial);
for ii=1 % For neatness, the pre-allocation variable list is hidden in this folded (dummy) FOR loop. 
    varList =   {
        
                 'cellID',      'string'; ...
                 'dataset',     'string'; ...     
                 'cellNo',      nan; ...
                 'trialNo',     nan(1,prms.maxNTrial); ...
                  
                'rMap',          scoreDum; ...
                 'SI',          scoreDum; ...
                 'intraStab',   scoreDum; ...

                     };
    varList = varList';
end
dummyRow = cell2table( varList(2,:) );
Res = dummyRow;
Res.Properties.VariableNames = varList(1,:);
Res.Properties.UserData = prms;

% Loop though results arrays from above and assign data to rows in table %
cellCount=1;
for ii=1:length(datasetNamesTemp)  %% ii = iterator for dataset
    
    for jj=1:size(siTemp{ii},2)  %% jj = iterator for cell
        
        Res(cellCount,:) = dummyRow;
        
        % Assign metadata for Cell %
        data = evalin('base',SD.selData{ii});
        Res.cellID(cellCount,1) = {[datasetNamesTemp{ii} ' t' num2str(data.trials(1).cells(jj).tet) 'c' num2str(data.trials(1).cells(jj).cellnum)]};
        Res.dataset(cellCount,1) = datasetNamesTemp(ii);
        Res.cellNo(cellCount,1) = jj;
        
        % Assign data by trial %
        for kk=1:size(siTemp{ii},1)  %% kk=iterator for trial
            
            Res.trialNo(cellCount,kk) = kk;
%             Res.rMap(cellCount,kk) = { squeeze( rMapTemp{ii}(kk,jj,:)) };
            Res.SI(cellCount,kk) = { squeeze( siTemp{ii}(kk,jj,:)) };
            Res.intraStab(cellCount,kk) = { squeeze( stabTemp{ii}(kk,jj,:)) };
            
        end
        
        % Bump the cell count %
        cellCount = cellCount + 1;
    end
    
end
            








    
