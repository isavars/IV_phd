function [shufT] = shuffleRMapsByCell_v2
% Spike shuffle rMaps. Output will be a data table with the same format as
% from 'getSpatData'. Fields for scores will have 1xprms.shufNSteps vectors
% of shuffled scores. Only 'HP' (baseline) trials will be shuffled.
% 03/03/2020 IV only famBox baseline trials to be shuffled.

%% params
prms.shufMinOffset = 20;       % In sec (minimum shift of spike train)
prms.shufNSteps = 1000;        % shuffles/cell (1000)
%map params
prms.mapType = 'boxcar'; %'boxcar'; 
prms.binSize = 10; %cam pix/spatial bin
% prms.adSmooth = 200;  
prms.smooth = 5;% size smoothing kernel (in bins)
prms.speedfilter = [2.5 400]; % lower/uppper limits cm/s

prms.dirBin = 6; % how many degress/bin
prms.dirSmooth = 5; %kernel size 
%path scaling
prms.minOccForEdge = 50; %n of pos samples (i.e. 1s)
prms.boxExtent = 250; %size map in cam pix
prms.mapSize = 25; % size of maps in bins
%trial selection
prms.TrialMode = 'fam_all'; %'HP_all'; %'HP_all'; 'HP_preProbe'
prms.maxNTrial = 5; %max n of trials of any data set - need to pre-determine
% if strcmp(prms.TrialMode,'HP_preProbe')
%     prms.maxNTrial = 4;
% end
% for border score
prms.binSizeCM = 2.5; %spatial bin size in cm
prms.rateThr = 0.3;       % Fraction of max rate to use as field threshold.
prms.sizeThr = 200;  
% type of gridness to use
prms.gridness = 'Sargolini'; % 'Wills'; 'Sargolini';
%save
prms.saveMaps = 1; %y/n
prms.rootDir = 'S:/DBIO_TKFC_SKGTIVA/';%'/home/laurenz/serverlink/lmuessig/!postDoc/!analysis/'; %'D:\'; %where to save maps? '/home/laurenz/serverlink/lmuessig/!postDoc/!analysis'

%% shuffling
% grab map params
if strcmp(prms.mapType,'adSmooth')
    rateMapParams_pups = rates_params('bin',prms.binSize,'adaptive_smooth',prms.adSmooth,'filt_speed',prms.speedfilter);
    rateMapParams_ad = rates_params('bin',prms.binSize,'adaptive_smooth',prms.adSmooth,'filt_speed',[prms.speedfilter(1)*2 prms.speedfilter(2)]);
elseif strcmp(prms.mapType,'boxcar')
    rateMapParamsSpat = rates_params('bin',prms.binSize,'smooth',prms.smooth,'filt_speed',prms.speedfilter);
    rateMapParamsDir = rates_params('space','dir','bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth);
    rateMapParamsDirPos = rates_params('space','dir','mode','pos','bin_dir',prms.dirBin,'smooth_dir',prms.dirSmooth);
end

SD = gss; %grab Scan datasets and load them into workspace
bigDataArray = cell(size(SD.selData,2),2);
for i=1:size(SD.selData,2)
    bigDataArray{i,1} = evalin('base',SD.selData{i});
    bigDataArray{i,2} = SD.selData{i};
end
%initiate the main scores 
[cellStr,dataSet,SI_spat,SI_dir,BS,gridness,RV_length] = deal({});

% Main loop for getting shuffled data scores %
for i=1:size(bigDataArray,1)
    data=bigDataArray{i,1};
    %get all trial types run with current data set
    env = cell(1,length(data.trials));
    for b=1:length(data.trials)
        env{b} = data.trials(b).user.environment;
    end
    trialInd = false(1,length(data.trials));
    
    %select appropriate trials
    if strcmp(prms.TrialMode,'fam_preProbe')
        nPreProbefam = find(~strcmp(env,'nov'),1,'first')-1;
        if nPreProbefam == 0
            nPreProbefam = length(data.trials);
        end
        trialInd(1:nPreProbefam) = true;
    elseif strcmp(prms.TrialMode,'fam_all')
        temp_idx = strcmp(env,'fam');
        trialInd(temp_idx) = true;
    end
    % get a column index so output will have sasme format as normal data
    % table
    colInd = 1:length(data.trials);
    colInd = colInd(trialInd);
    data.trials = data.trials(trialInd); % only select relevant trials
    
    %scale path
    [data_scaled,~] = lm_meanRateMap_findEdgesAndScalePath_CA1(data,prms.minOccForEdge,prms.boxExtent);
    %important for path scaling
    for c=1:length(data_scaled.trials)
        data_scaled.trials(c).window_x = 512;
        data_scaled.trials(c).window_y = 512;
    end
    
    % Unfortunately, to take account of different trial lengths, we have to split the analyis by trial %   
    [rMapTempDS,dMapTempDS] = deal(cell( length(data_scaled.trials), length(data_scaled.trials(1).cells), prms.shufNSteps));  
    [pMapTempDS,dPosMapsTempDS] = deal(cell( length(data_scaled.trials), 1, prms.shufNSteps));% Pre-allocate temp arrays for this dataset.
    [BSTemp, SI_spatTemp, SI_dirTemp, gridnessTemp, rvTemp] = deal(cell(  length(data_scaled.trials(1).cells),prms.maxNTrial));
    
    genTrialDur = min([data_scaled.trials(:).dur]); %use minimum duration of all trials in set as max 
    offsets = linspace( prms.shufMinOffset, genTrialDur - prms.shufMinOffset, prms.shufNSteps); %create linear spike train shifts
    
    % Also need to loop through *Each indivdual offset* %
    for k=1:length(offsets) %k = iterator for offset
        dataShuf = spk_randomise(data_scaled,'fixedWrap',offsets(k)); %shuffle data
        % Make maps 
        [rMaps, posMaps] = rates_main(dataShuf,rateMapParamsSpat);
         dMaps = rates_main(dataShuf,rateMapParamsDir);
         dPosMaps = rates_main(dataShuf,rateMapParamsDirPos);

        %store maps for later
        rMapTempDS(:,:,k) = rMaps;
        pMapTempDS(:,:,k) = posMaps; 
        dMapTempDS(:,:,k) = dMaps;
        dPosMapsTempDS(:,:,k) = dPosMaps;
        
        for t = 1:length(dataShuf.trials)
            % Get the scores for each cell
            for m=1:length(dataShuf.trials(t).cells)
%                 BSTemp{m,colInd(t)}(end+1) = lm_borderScoreForScaledRateMaps(rMaps{t,m}(1:prms.mapSize,1:prms.mapSize),prms.binSizeCM^2,'rateThr',prms.rateThr,'sizeThr',prms.sizeThr); %border score
                SI_spatTemp{m,colInd(t)}(end+1) = map_skaggsinfo(rMaps{t,m}(1:prms.mapSize,1:prms.mapSize),posMaps{t,1}(1:25,1:25) ); %spatial info
%                 SI_dirTemp{m,colInd(t)}(end+1) = map_skaggsinfo(dMaps{t,m},dPosMaps{t,1} ); %spatial info - dir
                %gridness
                sAC = map_crosscorr(rMaps{t,m}(1:prms.mapSize,1:prms.mapSize),rMaps{t,m}(1:prms.mapSize,1:prms.mapSize));
%                 if strcmp(prms.gridness,'Wills')
%                     gridnessTemp{m,colInd(t)}(end+1) = map_gridprops(sAC,'peakMode','point','corrThr',0,'radius','est'); 
%                 elseif strcmp(prms.gridness,'Sargolini')
% %                     binsize = (1/data_scaled.trials(t).ppm) * prms.binSize * 100;
%                     minPeakArea = 225/(prms.binSizeCM^2);
%                     gridnessTemp{m,colInd(t)}(end+1) = map_gridprops(sAC,'peakMode','area','areaThr',minPeakArea,'corrThr',0.1,'radius','fieldExtent'); %gridness
%                 end
                
                rvTemp{m,colInd(t)}(end+1) = dir_rayleighvector(dMaps{t,m}); %length of raleigh vector
                
                if k==1 && t==1
                    cellStr = cat(1,cellStr,{[bigDataArray{i,2} ' t' num2str(data.trials(1).cells(m).tet) 'c' num2str(data.trials(1).cells(m).cellnum)]}); %cell ID string
                end
            end
        end
    end
    % Transfer the temp arrays for this dataset to the overall cell array %
    dataSet = cat(1,dataSet,repmat(bigDataArray(i,2),m,1));
    SI_spat = cat(1,SI_spat,SI_spatTemp);
    SI_dir = cat(1,SI_dir,SI_dirTemp);
    RV_length = cat(1,RV_length,rvTemp);
    gridness = cat(1,gridness,gridnessTemp);
    BS = cat(1,BS,BSTemp);
    %save maps - good idea in case one needs to go back to them later
    if prms.saveMaps
        
        if ~isdir(prms.rootDir)
            mkdir(prms.rootDir);
        end
        
        rMapsCurrData = rMapTempDS;
        save([prms.rootDir bigDataArray{i,2} '_rMaps.mat'], 'rMapsCurrData');
        %check for file size problems
        [msgstr, ~] = lastwarn; %fetch warnings
        if ~isempty(msgstr)
            save([prms.rootDir  bigDataArray{i,2} '_rMaps.mat'],'-v7.3', 'rMapsCurrData');
        end
        
        dMapsCurrData = dMapTempDS;
        save([prms.rootDir  bigDataArray{i,2} '_dMaps.mat'], 'dMapsCurrData');
        %check for file size problems
        [msgstr, ~] = lastwarn; %fetch warnings
        if ~isempty(msgstr)
            save([prms.rootDir  bigDataArray{i,2} '_dMaps.mat'],'-v7.3', 'dMapsCurrData');
        end
        
        pMapsCurrData = pMapTempDS;
        save([prms.rootDir  bigDataArray{i,2} '_posMaps.mat'], 'pMapsCurrData');
        [msgstr, ~] = lastwarn; %fetch warnings
        if ~isempty(msgstr)
            save([prms.rootDir  bigDataArray{i,2} '_posMaps.mat'],'-v7.3', 'pMapsCurrData');
        end   
        
        dPosMapsCurrData = dPosMapsTempDS;
        save([prms.rootDir  bigDataArray{i,2} '_posMapsDir.mat'], 'dPosMapsCurrData');
        [msgstr, ~] = lastwarn; %fetch warnings
        if ~isempty(msgstr)
            save([prms.rootDir  bigDataArray{i,2} '_posMapsDir.mat'],'-v7.3', 'dPosMapsCurrData');
        end  
    end
    fprintf('Finished shuffling dataset %s\n',bigDataArray{i,2});
end
%
%% output
varList = {'dataset','cellID','SI_spatShuff','rvLengthShuff','SI_dirShuff','BorderScoreShuff','gridnessShuff'};
shufT = table(dataSet,cellStr,SI_spat,RV_length,SI_dir,BS,gridness,'VariableNames',varList);
shufT.Properties.UserData = prms;

save([prms.rootDir 'r889_FullDataSet_shufT'],'shufT');

end
            








    
