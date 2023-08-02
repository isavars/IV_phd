function [Res] = shuffledRVByCell(varargin)
% Spike shuffled thresholds for directionality for single rate maps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prms.shufMinOffset = 20;       % In sec
prms.shufNSteps = 1000;        %
prms.dirBin = 6;               %
prms.dirSmooth = 5;            %
prms.timeSliceWindow = 60;    % In sec
prms.maxNTrial = 10;
if ~isempty(varargin)
    for ii=1:(length(varargin)-1)
        prms.(varargin{ii}) = varargin{ii+1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SD=gss;
rateMapParams = rates_params('bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth,'space','dir');
posMapParams = rates_params('bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth,'space','dir','mode','pos');

for ii=1:length(SD.selData)
    bigDataArray{ii,1} = evalin('base',SD.selData{ii});
    bigDataArray{ii,2} = SD.selData{ii};
end

[rvTemp,siTemp,rvTempTS,datasetNamesTemp] = deal(cell(1,length(SD.selData)));   % Pre-allocate temp cell arrays for results (cannot combine into structure until after PARFOR).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop for getting shuffled data scores %
tic
parfor ii=1:length(SD.selData)
    data=bigDataArray{ii,1};
    
    % Unfortunately, to take account of different trial lengths, we have to spilt the analyis by trial %   
    [rvTempDS,siTempDS,stabTempDS,rvTempTSDS,stabTempTSDS] = deal(nan( length(data.trials), length(data.trials(1).cells), prms.shufNSteps));   % Pre-allocate temp arrays for this dataset.
    
    for jj=1:length(data.trials)  %% jj = iterator for trial (but note that the first action here is to reduce the data to one trial only, so after that jj only used for output assignment).
        trData=data;
        trData.trials = trData.trials(jj);
        offsets = linspace( prms.shufMinOffset, trData.trials.dur-prms.shufMinOffset, prms.shufNSteps);
        
        % Also need to loop through *Each indivdual offset* %
        for kk=1:length(offsets) %% kk = iterator for offset
            
            dataShuf = spk_randomise(trData,'fixedWrap',offsets(kk));
            
            % Make maps (just for 1 trial, one offset %
            rateMaps = rates_main(dataShuf,rateMapParams);
            posMaps = rates_main(dataShuf,posMapParams);
            tempParams = rateMapParams;   tempParams.filt_1stHalf = 1;
            maps1stHalf = rates_main(dataShuf,tempParams);
            tempParams.filt_1stHalf = 0;   tempParams.filt_2ndHalf = 1;
            maps2ndHalf = rates_main(dataShuf,tempParams);
            
            % Time slice maps %
            a = prms.timeSliceWindow;
            nMaps = floor(data.trials(jj).dur / a);
            timeLims = [0:a:((nMaps*a)-a); a:a:(nMaps*a)]';
            timeSliceMaps = cell(1,length(trData.trials.cells),nMaps);
            for mm=1:nMaps
                temp = rates_main(dataShuf, rates_params('space','dir','bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth,'filt_time',timeLims(mm,:)));
                timeSliceMaps(1,:,mm) = temp(1,:);
            end
            timeSliceRV = nanmean(      cellfun(@dir_rayleighvector,timeSliceMaps),     3  );
            timeSliceCorr = nanmean(       cellfun(@map_spatialcorr,timeSliceMaps(1,:,1:end-1),timeSliceMaps(1,:,2:end)),     3  );
            
            % Get the scores for each cell
            for mm=1:length(trData.trials.cells)
                dMap = rateMaps{1,mm};
                rvTempDS(jj,mm,kk) = dir_rayleighvector(dMap);
                siTempDS(jj,mm,kk) = map_skaggsinfo(dMap,posMaps{1});
                stabTempDS(jj,mm,kk) = map_spatialcorr( maps1stHalf{1,mm}, maps2ndHalf{1,mm} );
                rvTempTSDS(jj,mm,kk) = timeSliceRV(1,mm);
                stabTempTSDS(jj,mm,kk) = timeSliceCorr(1,mm);
            end

        end
        
        % Transfer the temp arrays for this dataset to the overall cell array %
        rvTemp{ii} = rvTempDS;
        siTemp{ii} = siTempDS;
        stabTemp{ii} = stabTempDS;
        rvTempTS{ii} = rvTempTSDS;
        stabTempTS{ii} = stabTempTSDS;
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

                 'RV',          scoreDum; ...
                 'SI',          scoreDum; ...
                 'intraStab',   scoreDum; ...
                 'timeSliceRV', scoreDum; ...
                 'timeSliceStab', scoreDum; ...

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
    
    for jj=1:size(rvTemp{ii},2)  %% jj = iterator for cell
        
        Res(cellCount,:) = dummyRow;
        
        % Assign metadata for Cell %
        data = evalin('base',SD.selData{ii});
        Res.cellID(cellCount,1) = {[datasetNamesTemp{ii} ' t' num2str(data.trials(1).cells(jj).tet) ' c' num2str(data.trials(1).cells(jj).cellnum)]};
        Res.dataset(cellCount,1) = datasetNamesTemp(ii);
        Res.cellNo(cellCount,1) = jj;
        
        % Assign data by trial %
        for kk=1:size(rvTemp{ii},1)  %% kk=iterator for trial
            
            Res.trialNo(cellCount,kk) = kk;
            Res.RV(cellCount,kk) = { squeeze( rvTemp{ii}(kk,jj,:)) };
            Res.SI(cellCount,kk) = { squeeze( siTemp{ii}(kk,jj,:)) };
            Res.intraStab(cellCount,kk) = { squeeze( stabTemp{ii}(kk,jj,:)) };
            Res.timeSliceRV(cellCount,kk) = { squeeze( rvTempTS{ii}(kk,jj,:)) };
            Res.timeSliceStab(cellCount,kk) = { squeeze( stabTempTS{ii}(kk,jj,:)) };
            
        end
        
        % Bump the cell count %
        cellCount = cellCount + 1;
    end
    
end
            








    
