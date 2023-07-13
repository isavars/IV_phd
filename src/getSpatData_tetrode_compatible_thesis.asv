function [ResT] = getSpatData_tetrode_compatible_thesis(varargin)
% Finds place and directional spatial firing properties for all cells in
% loaded scan data.TO DO - 1) trimm table down to only what is needed for
% thesis analysis, 2) add variable trial lengths and the condition for
% trial 6 needs to be for trials named 'sleep' instead. 

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

% Scale data to fit into standardised square.
prms.scaleToFitSquare = 0;
prms.minOccForEdge    = 50; %in pos samples (i.e. 1s)
prms.boxExtent        = 250; %n of cam pix for 62.5 cm box

prms.binSizeCM = 2.5; %Spatial bin size in cm - for border score
prms.rateThr = 0.3;   %Fraction of max rate to use as field threshold.
prms.sizeThr = 200;  
% type of gridness to use
prms.gridness = 'Sargolini'; % 'Wills'; 'Sargolini';

prms.envToInclude = {'all'};
prms.maxNTrial    = 6;     % The most trials in any dataset. Needed for pre-allocation. IV changed to 6 . 
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
                 'trialName',   categorical(scoreDum); ... 
                 'trialDur',    scoreDum; ...
                 'env',         categorical(scoreDum); ...

                 'meanRate',    scoreDum; ...
                 'nSpks',       scoreDum; ...
                 'SpkTs',       cell(size(scoreDum)); ...
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
                 'waveforms',    cell(size(scoreDum)); ... %TW added
                 'wf_means',  cell(size(scoreDum)); ... %IV added 
                 

                 
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

    trialsToUse = 1:length(data.trials);
    envList = {};
   
    for jj = 1:length(data.trials)
        envList{jj} = data.trials(jj).user.environment;
    end
       
    % scale path (if required)
    if prms.scaleToFitSquare
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
          if isempty(ResT.env(cellCount,kk)) || ~strcmp(string(ResT.env(cellCount,kk)), 'sleep') %|| ResT.env(cellCount,kk) ~= envList(6)  
            if prms.scaleToFitSquare
                [RMWinR,RMWinC] = deal( 1:25 );
            else
                [RMWinR,RMWinC] = deal( 1:45 );%size (rMap); %if size is different code won't run error says different size and position input. runs if new indices are assigned everywhere rMap is used. 
            end
            
            dMap=dMaps{trialsToUse(kk),jj}; % dir map
            rMap=rMaps{trialsToUse(kk),jj}(RMWinR,RMWinC); % rate map
            
            %adding trial names and trial duration - needs to be in this loop
            ResT.trialName(cellCount,kk) = {data.trials(kk).trialname};
            ResT.trialDur(cellCount,kk) = data.trials(kk).dur;
        
%             if nansum(dMap(:))==0;   continue;   end            
            %%% Store the Rate map for later reference %%%
            ResT.dMap{cellCount,kk} = dMap;
            ResT.dMap1stHalf{cellCount,kk} = dMaps1stHalf{trialsToUse(kk),jj};
            ResT.dMap2ndHalf{cellCount,kk} = dMaps2ndHalf{trialsToUse(kk),jj};
            ResT.rMap{cellCount,kk} = rMap;
            ResT.rMap1stHalf{cellCount,kk} = rMaps1stHalf{trialsToUse(kk),jj}(RMWinR,RMWinC);
            ResT.rMap2ndHalf{cellCount,kk} = rMaps2ndHalf{trialsToUse(kk),jj}(RMWinR,RMWinC);
            
            %%% spatial Stats %%%    
            ResT.SI_spat(cellCount,kk) = map_skaggsinfo(rMap,rPosMaps{trialsToUse(kk),1}(RMWinR,RMWinC)); % spatial info - spat
            ResT.intraStabSpat(cellCount,kk) = map_spatialcorr( rMaps1stHalf{trialsToUse(kk),jj}(RMWinR,RMWinC) , rMaps2ndHalf{trialsToUse(kk),jj}(RMWinR,RMWinC)); %intra trial stability
            ResT.borderScore(cellCount,kk) = lm_borderScoreForScaledRateMaps(rMap,prms.binSizeCM^2,'rateThr',prms.rateThr,'sizeThr',prms.sizeThr); %border score   
            
            %%% Inter-trial comparisons %%%
            if kk+1<=length(trialsToUse) && trialsToUse(kk)+1==trialsToUse(kk+1)  % Only correlate temporally adjacent trials
                ResT.interStabDir(cellCount,kk) = map_spatialcorr(dMap,dMaps{trialsToUse(kk+1),jj}); %across trial stability - dir
                ResT.interStabSpat(cellCount,kk) = map_spatialcorr(rMap,rMaps{trialsToUse(kk+1),jj}(RMWinR,RMWinC)); % across trial stability - spat
            end  
            %%% Basic measures of firing %%%
            ResT.peakRate(cellCount,kk) = nanmax(rMap(:));  % Peak Rate rate map
            ResT.meanRate(cellCount,kk) = length( data.trials( trialsToUse(kk) ).cells(jj).st ) / data.trials( trialsToUse(kk) ).dur; % mean rate 
            ResT.nSpks(cellCount,kk) = length(data.trials( trialsToUse(kk) ).cells(jj).st ); % number of spikes
            ResT.SpkTs{cellCount,kk} = data.trials( trialsToUse(kk) ).cells(jj).st;
            ResT.burstIndex(cellCount,kk) = (sum(diff(data.trials( trialsToUse(kk) ).cells(jj).st) <= 0.009))/(length(diff(data.trials( trialsToUse(kk) ).cells(jj).st))); % burst index is the proportion of all inter spike intervals that were <= 6ms - IV             
            
            % Get the mean waveform (just channel with highest amplitude).
            if ~isempty(data_scaled.trials(trialsToUse(kk)).cells(jj).wf_means) % This contains the mean waveform for each cluster, format (1:50,1:4), time x tetrode channel.
                wf = (data_scaled.trials(trialsToUse(kk)).cells(jj).wf_means);
            else
                wf = NaN;
            end
            wf           = (wf./128) .* (data_scaled.trials(trialsToUse(kk)).cells(jj).scalemax);  % Convert to true voltage (units in scan are -127:128 digitisation levels).
            wfMins       = nanmin( wf, [], 1 );
            wfMaxs       = nanmax( wf, [], 1 );
            wfAmps       = wfMaxs - wfMins;
            [~,maxAmpCh] = nanmax( wfAmps );
            ResT.waveforms{cellCount,kk} = wf( :, maxAmpCh );
            ResT.wf_means{cellCount,kk} = (data_scaled.trials(trialsToUse(kk)).cells(jj).wf_means);
            
        %Filter sleep trial for State - trying out with no sleep filter to
        %make wfs
        elseif strcmp(string(ResT.env(cellCount,kk)), 'sleep') 

            %sleep option 1 
%            [bestEEG, peakTheta, peakDelta] = select_eeg; 
%            stateInds =  getBrainStateHardThr (data.trials(trialsToUse(kk)).speed, bestEEG, 'thetaFr', peakTheta, 'deltaFr', peakDelta);
%            SWS = stateInds.sws;
%            indStart = diff([0,SWS,0]) == 1;  indStart = find(indStart == 1); indStart = indStart*0.8 - 0.4;
%            indStop = diff([0,SWS,0]) == -1; indStop = find(indStop == 1); indStop = indStop*0.8 + 0.4; 
%            epochInds = [indStart; indStop];
%             spikeTimes = data.trials(trialsToUse(kk)).cells(jj).st; 
%            newSpikeTimes = [];
%            newTrialDuration = sum(diff(epochInds));
%            for itEps=1:size(epochInds, 2)
%                indEps = spikeTimes >= epochInds(1,itEps) & spikeTimes <= epochInds(2,itEps);
%                newSpikeTimes = [newSpikeTimes; spikeTimes(indEps,:)];   
%            end
%             %%% Basic measures of firing %%%
%             ResT.peakRate(cellCount,kk) = nanmax(rMap(:));  % Peak Rate rate map
%             ResT.meanRate(cellCount,kk) = length(newSpikeTimes) / newTrialDuration; % mean rate 
%             ResT.nSpks(cellCount,kk) = length(newSpikeTimes ); % number of spikes
%             ResT.SpkTs{cellCount,kk} = newSpikeTimes;
%             ResT.burstIndex(cellCount,kk) = (sum(diff(newSpikeTimes) <= 0.009))/(length(diff(newSpikeTimes))); % burst index (Knierim criteria) is the proportion of all inter spike intervals that were <= 6ms - IV 
%             % Sleep adjusted trial duration 
%             ResT.trialDur(cellCount,kk) = newTrialDuration;
            
            %sleep option 2- for making spatData to be used in wf making 
            %adding trial names and trial duration - not adjusted for sleep
            spikeTimes = data.trials(trialsToUse(kk)).cells(jj).st; 
            ResT.trialName(cellCount,kk) = {data.trials(kk).trialname};
            ResT.trialDur(cellCount,kk) = data.trials(kk).dur;
            ResT.trialName(cellCount,kk) = {data.trials(kk).trialname};
            %%% Basic measures of firing %%% - not adjusted for sleep 
            ResT.peakRate(cellCount,kk) = nanmax(rMap(:));  % Peak Rate rate map
            ResT.meanRate(cellCount,kk) = length(spikeTimes) / data.trials(kk).dur; % mean rate 
            ResT.nSpks(cellCount,kk) = length(spikeTimes ); % number of spikes
            ResT.SpkTs{cellCount,kk} = spikeTimes;
            ResT.burstIndex(cellCount,kk) = (sum(diff(spikeTimes) <= 0.009))/(length(diff(spikeTimes)));

            % Get the mean waveform (just channel with highest amplitude).
            if ~isempty(data_scaled.trials(trialsToUse(kk)).cells(jj).wf_means) % This contains the mean waveform for each cluster, format (1:50,1:4), time x tetrode channel.
                wf = (data_scaled.trials(trialsToUse(kk)).cells(jj).wf_means);
            else
                wf = NaN;
            end
            wf           = (wf./128) .* (data_scaled.trials(trialsToUse(kk)).cells(jj).scalemax);  % Convert to true voltage (units in scan are -127:128 digitisation levels).
            wfMins       = nanmin( wf, [], 1 );
            wfMaxs       = nanmax( wf, [], 1 );
            wfAmps       = wfMaxs - wfMins;
            [~,maxAmpCh] = nanmax( wfAmps );
            ResT.waveforms{cellCount,kk} = wf( :, maxAmpCh );
            ResT.wf_means{cellCount,kk} = (data_scaled.trials(trialsToUse(kk)).cells(jj).wf_means);

          end                         
        end
        
        % Bump up the cell count %
        cellCount = cellCount + 1;
    end
end
close(hWait);
% successSig;
end

function [bestEEG, peakTheta, peakDelta] = select_eeg(varargin)
%SELECT EEG compares eegs by making power spectrums from a scan file and selects the best one 
%   Takes in data from scan file and outputs bestEEG
%   TO DO: 
%       1. ADD exclude eegs with an snr of under 3 
%       2. speed filter the theta - use speed and position data -
%       each speed point for 5 eeg points - make logical array of indexes
%       (filtering for speeds greater than 2) %01/06/23 Iv - is this fixed?
%       

SD=gss;%grab all scan datasets 
for itSD =1:length(SD.selData)
    data = evalin('base', SD.selData{itSD}); %load current data set into work space
    for  itTR = 1:length(data.trials) % make indexes for baseline and sleep 
        if strcmp(data.trials(itTR).user.environment, 'fam')
            baselineTrial = itTR;
        elseif strcmp(data.trials(itTR).user.environment, 'sleep')
            sleepTrial = itTR;
        end
    end
    
% get best eeg and delta frequencies from sleep trial: 

    %get eegs for sleep trial and baseline trial from current folder. 
    sleeptrialname = [data.load_selection.pathName data.trials(sleepTrial).trialname '.eeg'];
    baselinetrialname = [data.load_selection.pathName data.trials(baselineTrial).trialname '.eeg'];
    eegs_S = [];
    eegs_B = [];
    for it_eegs = 1:32
        if it_eegs == 1
            eeg_struct_s  = load_eeg(sleeptrialname);
            eeg_struct_b  = load_eeg(baselinetrialname);
        else
            eeg_struct_s = load_eeg(strcat(sleeptrialname,num2str(it_eegs)));
            eeg_struct_b = load_eeg(strcat(baselinetrialname,num2str(it_eegs)));
        end
        eegs_S = [eeg_struct_s.eeg;eegs_S];
        eegs_B = [eeg_struct_b.eeg;eegs_B];
    end
    eegs_S = reshape(eegs_S,[],32);
    eegs_B = reshape(eegs_B,[],32);

    SNR = [];    
    for itEEG = 1:32 % hc maybe find something that does channel numbers 
        [~, ~, snr] = eeg_powerspec(eegs_S(:,itEEG),data.trials(sleepTrial).eeg(1).sample_rate); 
        SNR(itEEG) = snr;
    end
    [~,bestSNR] = max(SNR); %this finds the index for the eeg with the best SNR 
    peakFreqD = eeg_powerspec(eegs_S(:,bestSNR),data.trials(sleepTrial).eeg(1).sample_rate,'thetaBand',[1.5 4],'hfCutOff',25);

% get peak theta from baseline trial:
    speedFilteredEEG = [];
    baselineEEG = eegs_B(:,bestSNR);
    filteredSpeed = data.trials(baselineTrial).speed >= 2;% this needs to account for the 5 to 1 eeg to speed 
    filteredSpeed = reshape(repmat(filteredSpeed,1,5).',1,[]);
    for itSp = 1: length(filteredSpeed)
         if filteredSpeed(itSp) == 1 
            speedFilteredEEG = [speedFilteredEEG; baselineEEG(itSp)];
         end
    end
    peakFreqT = eeg_powerspec(speedFilteredEEG,data.trials(baselineTrial).eeg(1).sample_rate,'thetaBand',[5 11],'hfCutOff',15);

% outputs 
    bestEEG = eegs_S(:,bestSNR);
    peakTheta = peakFreqT;
    peakDelta = peakFreqD;
end
end