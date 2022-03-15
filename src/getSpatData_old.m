 function [ResT] = getSpatData(varargin)
% Place and directional spatial firing properties %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis parameters %
% maps - dir
prms.dirBin       = 6; % how many degress/bin
prms.dirSmooth    = 5; % kernel size
% maps - place
prms.placeBin     = 10; %cam pix/spatial bin
prms.placeSmooth  = 5;  % size smoothing kernel (in bins)
prms.adapt_smooth = 200; % 
prms.crop         = [450 450];
prms.speedFilt    = [2.5 400]; % lower/uppper limits cm/s
prms.minOccForEdge = 50; %in pos samples (i.e. 1s)
prms.boxExtent = 250; %n of cam pix for 62.5 cm box

prms.binSizeCM = 2.5; %Spatial bin size in cm - for border score
prms.rateThr = 0.3;   %Fraction of max rate to use as field threshold.
prms.sizeThr = 200;  
% type of gridness to use
prms.gridness = 'Sargolini'; % 'Wills'; 'Sargolini';

prms.envToInclude = {'all'};
prms.maxNTrial    = 5;     % The most trials in any dataset. Needed for pre-allocation. IV changed to 5. 
if ~isempty(varargin)
    for ii=1:(length(varargin)-1)
        prms.(varargin{ii}) = varargin{ii+1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up results table %
scoreDum = nan(1,prms.maxNTrial);
for ii=1 % For neatness, the pre-allocation variable list is hidden in this folded (dummy) FOR loop. 
    varList =   {
                 'cellID',      'string'; ...
                 'animal',      'string'; ...   
                 'dataset',     'string'; ...     
                 'cellNo',      nan; ...
                 'trialNo',     scoreDum; ...
                 'env',         categorical(scoreDum); ...

                 'meanRate',    scoreDum; ...
                 'nSpks',       scoreDum; ...
                 'peakRate',    scoreDum; ...
                 'burstIndex',  scoreDum; ... %IV added 

                 'pl',          scoreDum; ...
                 'meanSpeed',   scoreDum; ...
                 'nBinVis',     scoreDum; ...

                 'RV_length',   scoreDum; ...
                 'RV_dir',      scoreDum; ...
                 'SI_dir',          scoreDum; ...
                 'interStabDir',   scoreDum; ...
                 'intraStabDir',   scoreDum; ...
                 
                 'SI_spat',          scoreDum; ...
                 'interStabSpat',   scoreDum; ...
                 'intraStabSpat',   scoreDum; ...
                 'borderScore',     scoreDum; ...
                 
                 'sAC',         cell(size(scoreDum)); ...
                 'gridness'     scoreDum; ...
                 'wavelength'     scoreDum; ...
                 'orientation'     scoreDum; ...
                 
                 'rMap',     cell(size(scoreDum)); ...
                 'rMap1stHalf', cell(size(scoreDum)); ... 
                 'rMap2ndHalf', cell(size(scoreDum)); ... 
                 'dMap',     cell(size(scoreDum)); ...
                 'dMap1stHalf', cell(size(scoreDum)); ... 
                 'dMap2ndHalf', cell(size(scoreDum)); ... 
                 'waveform',    cell(size(scoreDum)); ...
                 
                     };
    varList = varList';
end
ResT = cell2table( varList(2,:) );
ResT.Properties.VariableNames = varList(1,:);
ResT.Properties.UserData = prms;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main analysis loop %
cellCount = 1;
SD=gss; %grab all Scan datasets
hWait=waitbar(0);
for ii=1:length(SD.selData)
    data = evalin('base', SD.selData{ii}); %load current data set into work space
    waitbar(ii/length(SD.selData),hWait,SD.selData{ii});
    
    %%% Get the list of trials to analyse %%%
    
%     trialsToUse  = [];
%     envList = {};
%     if strcmp( prms.envToInclude{1}, 'all' )
%         trialsToUse = 1:length( data.trials );
%         for jj=1:length(data.trials) 
%             envList{jj}=data.trials(jj).user.environment; %#ok<AGROW>
%         end
%     else
%         for jj=1:length(data.trials)
%             if any(strcmp(prms.envToInclude, data.trials(jj).user.environment))
%                 trialsToUse(end+1) = jj; %#ok<AGROW>
%                 envList(end+1) = {data.trials(jj).user.environment}; %#ok<AGROW>
%             end
%         end
%     end

    trialsToUse = [];
    envList = {};
    trialOrder = {'fam', 'fam', 'nov', 'fam', 'sleep'};
    if strcmp(prms.envToInclude{1}, 'all' )
        for jj = 1:length(data.trials)
            envList{jj} = data.trials(jj).user.environment;
        end
        for jj = 1:length(data.trials)
            if ~any(strcmp(envList(jj), trialOrder(jj)))
            NewEnvList = cell(1,5);
                if strcmp(data.trials(2).user.environment, trialOrder(3))
                    for kk= length(data.trials):-1:1
                        data.trials (kk+1) = data.trials(kk);
                    end
                    NewEnvList{1} = ''; 
                    NewEnvList(2:end) = envList;
        %                     data.trials(1) = struct(1,length(fieldnames(data.trials)))%insert empty contents for data.trials(1)
                elseif strcmp(data.trials(4).user.environment,trialOrder(5))
                    data.trials(5) = data.trials(4);
        %                     emptyRow = cell(1,length(fieldnames(data.trials)))
        %                     data.trials(4) = cell2struct(emptyRow) %%insert empty contents for data.trials(4)
                    NewEnvList{4} = '';                   
                    NewEnvList{5} = envList{end};
                    NewEnvList(1:3) = envList(1:3);
                else
                    strcmp(data.trials(1).user.environment,trialOrder(5));
                    data.trials(5)= data.trials(1);
                    data.trials(1)= data.trials(2);
                    NewEnvList{5} = envList{1};
                    NewEnvList{1} = '';
                    NewEnvList(2:4) = envList(2:4);   
                end
            envList = NewEnvList;
            end
        end
        trialsToUse = 1:length(data.trials);
    end  
    
    % scale path (if required)
    if 0
        [data_scaled,~] = lm_meanRateMap_findEdgesAndScalePath_v2(data,prms.minOccForEdge,prms.boxExtent);
    else
       data_scaled = data; 
    end
    
    %important for path scaling
    for c=1:length(data_scaled.trials)
        data_scaled.trials(c).window_x = 512;
        data_scaled.trials(c).window_y = 512;
    end
    %%% Make maps %%%
    %dir maps
    dMaps = rates_main(data_scaled, rates_params('space','dir','bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth,'trial_index',trialsToUse));
    dMaps1stHalf = rates_main(data_scaled, rates_params('space','dir','bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth,'filt_1stHalf',1,'trial_index',trialsToUse));
    dMaps2ndHalf = rates_main(data_scaled, rates_params('space','dir','bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth,'filt_2ndHalf',1,'trial_index',trialsToUse));
    dPosMaps = rates_main(data_scaled, rates_params('space','dir','mode','pos','bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth,'trial_index',trialsToUse));
    %place maps
    [rMaps,rPosMaps] = rates_main(data_scaled,rates_params('bin',prms.placeBin,'smooth',prms.placeSmooth,'filt_speed',prms.speedFilt,'trial_index',trialsToUse,'adaptive_smooth',prms.adapt_smooth,'crop',prms.crop));
    rMaps1stHalf = rates_main(data_scaled,rates_params('bin',prms.placeBin,'smooth',prms.placeSmooth,'filt_speed',prms.speedFilt,'filt_1stHalf',1,'trial_index',trialsToUse,'adaptive_smooth',prms.adapt_smooth,'crop',prms.crop));
    rMaps2ndHalf = rates_main(data_scaled,rates_params('bin',prms.placeBin,'smooth',prms.placeSmooth,'filt_speed',prms.speedFilt,'filt_2ndHalf',1,'trial_index',trialsToUse,'adaptive_smooth',prms.adapt_smooth,'crop',prms.crop));
    
    %%% Analyse Path/Position/Behaviour data. This takes the form of one datapoint per trial, 
    %%% so we will calculate scores here, and then assign to every cell in the main loop below.
    pathScoresForDS = nan(3,prms.maxNTrial);  % A note on format: we store the path-type scores in a 4xnTrial cell array, which will get assigned into the Res table in the appropriate row for the cell, in the main loop below.
    for jj=1:length(trialsToUse)
        dist = sqrt( diff(double(data_scaled.trials(trialsToUse(jj)).x)).^2 + diff(double(data_scaled.trials(trialsToUse(jj)).y)).^2 );
        pathScoresForDS(1,jj) = (sum(dist)) / data_scaled.trials(trialsToUse(jj)).ppm;   % Row 1 = path length
        pathScoresForDS(2,jj) = mean(data_scaled.trials(trialsToUse(jj)).speed);         % 2 = mean speed
        pathScoresForDS(3,jj) = sum(sum(~isnan(rPosMaps{trialsToUse(jj),1})));       % 3 = N Position Bins Visited
    end 

    %%% Main data collection loop %%%
    for jj=1:length(data_scaled.trials(1).cells)   %% In main loop, jj = iterator for cell %%
        
        % Pre-assign the row for this cell %
        ResT(cellCount,:) = cell2table( varList(2,:) );
        
        %%% General metadata for the cell %%%
        ResT.cellID(cellCount,1) = {[SD.selData{ii} ' t' num2str(data.trials(1).cells(jj).tet) 'c' num2str(data.trials(1).cells(jj).cellnum)]};%cell ID string
        ResT.animal(cellCount,1) = cellstr(data.user.ID); %animal ID
        ResT.dataset(cellCount,1) = SD.selData(ii); %dataset name
        ResT.cellNo(cellCount,1) = jj; %cell number

        %%% Path and behaviour: assign for all cells here %%%
        ResT.pl(cellCount,:) = pathScoresForDS(1,:); %pathlengths
        ResT.meanSpeed(cellCount,:) = pathScoresForDS(2,:); % mean speed
        ResT.nBinVis(cellCount,:) = pathScoresForDS(3,:); %n of visited bins
        % Also assign some trial specific metadata %
        ResT.trialNo(cellCount,1:length(trialsToUse)) = trialsToUse;
        ResT.env(cellCount,1:length(envList)) = envList; %environment
        
        
        
        for kk=1:length(trialsToUse)             %% In main loop, kk = iterator for trial %%
            
            dMap=dMaps{trialsToUse(kk),jj}; % dir map
            rMap=rMaps{trialsToUse(kk),jj}(1:25,1:25); % rate map  
%             if nansum(dMap(:))==0;   continue;   end            
            %%% Store the Rate map for later reference %%%
            ResT.dMap{cellCount,kk} = dMap;
            ResT.dMap1stHalf{cellCount,kk} = dMaps1stHalf{trialsToUse(kk),jj};
            ResT.dMap2ndHalf{cellCount,kk} = dMaps2ndHalf{trialsToUse(kk),jj};
            ResT.rMap{cellCount,kk} = rMap;
            ResT.rMap1stHalf{cellCount,kk} = rMaps1stHalf{trialsToUse(kk),jj}(1:25,1:25);
            ResT.rMap2ndHalf{cellCount,kk} = rMaps2ndHalf{trialsToUse(kk),jj}(1:25,1:25);
            
            %%% Basic measures of firing %%%
            ResT.peakRate(cellCount,kk) = nanmax(rMap(:));  % Peak Rate rate map
            ResT.meanRate(cellCount,kk) = length( data.trials( trialsToUse(kk) ).cells(jj).st ) / data.trials( trialsToUse(kk) ).dur; % mean rate 
            ResT.nSpks(cellCount,kk) = length( data.trials( trialsToUse(kk) ).cells(jj).st ); % number of spikes
            ResT.burstIndex(cellCount,kk) = (sum(diff(data.trials( trialsToUse(kk) ).cells(jj).st) <= 0.009))/(length(diff(data.trials( trialsToUse(kk) ).cells(jj).st))); 
            % burst index is the proportion of all inter spike intervals that were <= 6ms - IV 
                       
            %%% Directional Stats %%%    
%             ResT.SI_dir(cellCount,kk) = map_skaggsinfo(dMap,dPosMaps{trialsToUse(kk),1}); % spatial info - dir
%             [ResT.RV_length(cellCount,kk), ResT.RV_dir(cellCount,kk)] = dir_rayleighvector(dMap); %raleigh vector props
%             ResT.intraStabDir(cellCount,kk) = map_spatialcorr( dMaps1stHalf{trialsToUse(kk),jj}, dMaps2ndHalf{trialsToUse(kk),jj} ); %intra trial stability
%             
            %%% spatial Stats %%%    
            ResT.SI_spat(cellCount,kk) = map_skaggsinfo(rMap,rPosMaps{trialsToUse(kk),1}(1:25,1:25) ); % spatial info - spat
            ResT.intraStabSpat(cellCount,kk) = map_spatialcorr( rMaps1stHalf{trialsToUse(kk),jj}(1:25,1:25) , rMaps2ndHalf{trialsToUse(kk),jj}(1:25,1:25) ); %intra trial stability
            ResT.borderScore(cellCount,kk) = lm_borderScoreForScaledRateMaps(rMap,prms.binSizeCM^2,'rateThr',prms.rateThr,'sizeThr',prms.sizeThr); %border score
 
            %%% grid Stats %%%
%             [ResT.sAC{cellCount,kk}] = map_crosscorr(rMap, rMap);
%             if strcmp(prms.gridness,'Wills')
%                 [ResT.gridness(cellCount,kk), Props] = map_gridprops(ResT.sAC{cellCount,kk},'peakMode','point','corrThr',0,'radius','est'); %gridness
%             elseif strcmp(prms.gridness,'Sargolini')
% %                 binsize = (1/data_scaled.trials(trialsToUse(kk)).ppm) * prms.placeBin * 100;
%                 minPeakArea = 225/(prms.binSizeCM^2);
%                 [ResT.gridness(cellCount,kk), Props] = map_gridprops(ResT.sAC{cellCount,kk},'peakMode','area','areaThr',minPeakArea,'corrThr',0.1,'radius','fieldExtent'); %gridness
%             end
%             ResT.wavelength(cellCount,kk) = Props.waveLength * prms.binSizeCM; %scale in cm
%             ResT.orientation(cellCount,kk) = Props.orientation; %orientation
%             
            %%% Inter-trial comparisons %%%
            if kk+1<=length(trialsToUse) && trialsToUse(kk)+1==trialsToUse(kk+1)  % Only correlate temporally adjacent trials
                ResT.interStabDir(cellCount,kk) = map_spatialcorr(dMap,dMaps{trialsToUse(kk+1),jj}); %across trial stability - dir
                ResT.interStabSpat(cellCount,kk) = map_spatialcorr(rMap,rMaps{trialsToUse(kk+1),jj}(1:25,1:25)); % across trial stability - spat
            end
            
            % Get the largest waveform.
%             figure; plot(
%             r889_191204_P18.trials(1).cells(2).wf_means(:,3) ); %this
%             prouces a waveform plot (don't use here but notes from what
%             Tom did 05/03/2020
% amps = max( r889_191204_P18.trials(1).cells(2).wf_means )  -  min( r889_191204_P18.trials(1).cells(2).wf_means )
% 
% amps =
% 
%   1�4 single row vector
% 
%   233.4265  233.3088  232.1324  233.2059
% 
% [~,maxCh] = max(amps)
% 
% maxCh =
% 
%      1

%more notes from the output of what tom did to get the max channel for
%waveform
            
        end
        
        % Bump the cell count %
        cellCount = cellCount + 1;
    end
end
close(hWait);
% successSig;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


