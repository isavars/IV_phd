function [Results,motionSensor] = logger_wrapper_all_files(varargin)
spikeThreshold = 60;
%% EDIT: File paths
% Ensure most recent masterfile path in CD before proceeding
% Home computer
trialInfo.masterfilePath = 'C:\Users\Tara\Tara Dropbox\Tara O''Driscoll\Logger project\MasterfileLogger.xlsx'; % Excel sheet with trial information
addpath('C:\Users\Tara\Tara Dropbox\Tara O''Driscoll\Logger project\Current logger code'); % Dropbox folder with code 

% Lab computer
% trialInfo.masterfilePath = 'C:\Users\Tara\Dropbox\Logger project\MasterfileLogger.xlsx';
% addpath('C:\Users\Tara\Dropbox\Logger project\Current logger code'); % Dropbox folder with code 

trialInfo.NDF_filepath = 'F:\LoggerProject\'; % Directory where NDF folders stored
trialInfo.DACQ_filepath = 'F:\DACQ files\'; % Directory where DACQ folders stored
trialInfo.BonsaiDir = 'F:\LoggerProject\BonsaiOutput'; % Bonsai output directory for pos data F:\LoggerProject\BonsaiOutput\Old name Bonsai files r861 r862
trialInfo.writeDir = 'F:\LoggerProject\TestingCodeBug';
% trialInfo.writeDir = 'F:\LoggerProject\Assess_Recording_Data';

%% Check for varargin
if strcmp(varargin{1},'Tint')
    process_pos = 0; % don't process position data if using Tint
elseif isempty(varargin{1})
    process_pos = 1;
end

if ~isempty(varargin{2})
    trialInfo.maDLC = 1;
end

%% Extract file names, trial info and pos data

[trialInfo,EventRecords] = readMasterFile(trialInfo);

if isempty(trialInfo.threshold)
    trialInfo.threshold = spikeThreshold;
end

fprintf('#################################################\n')
formatSpec = "Will search for logger files in folder path: %s\n";
fprintf(formatSpec,trialInfo.LogFilePath{1});

% Retrieve trial filenames 
[fileNames,trialInfo,EventInfo] = get_trial_files(trialInfo,EventRecords);

% Load DACQ trial settings from chosen set header
[DACQinfo,settxt] = getDACQinfo(trialInfo);

% Collect and post-process position data
[posInfo] = processTrackingData(trialInfo,process_pos,EventInfo);
EventInfo.numTouches = numel(EventInfo.TimeEventStart);
save(fullfile(trialInfo.writeDir,strcat(trialInfo.finalTrialName,'_EventInfo')),'EventInfo')

cd(trialInfo.LogFilePath{1}); % Look in folder with NLE/CSV file for logger files 
fprintf('#################################################\nLogger filenames retrieved....\n#################################################\n')

%% Generate Tint .set file
fprintf('Generating Set file....\n#################################################\n')
TintSetfile_header(trialInfo,settxt,posInfo);

%% Create Pos file 
fprintf('Generating Pos file....\n#################################################\n')
generatePos(trialInfo,posInfo,DACQinfo)

%% Preallocate for spike data 

new_tets = cell(4,numel(fileNames));
num_spikes = cell(4,numel(fileNames));
num_spikes(:,:) = {0};
all_LFP = cell(numel(fileNames),1);

if trialInfo.fileFormat == 2
    BlockDataAll = cell(4,numel(fileNames));
end

     
%% Parse through all files
jj = 1;
tic
while jj <= numel(fileNames)
    fprintf(strcat("Processing file ",num2str(jj)," of ",num2str(numel(fileNames)),'\n'))  
    if jj == 1
        time_start = 0; % beginning of trial
        time_end = ((2^19) - EventInfo.LoggerSamplesToRemove) * 32 - 2; % accounts for delay in Bonsai start
    else
        EventInfo.LoggerSamplesToRemove = []; % Don't remove samples
        time_start = time_end + 2;
        time_end = time_end + ((2^19) * 32 - 1); % Update time with each file
    end
    
    % Looks at current tetrode data with respect to other tetrodes (for correct
    % ampsamp)
    time_diff = 2; % 2 microseconds between each sample
    timestamps = time_start:time_diff:time_end;
         
      if trialInfo.fileFormat == 1 % flat file format
        [microvolts,EventInfo] = convert_to_volts(fileNames{jj}, EventInfo); % convert file to voltages  
      elseif trialInfo.fileFormat == 2 % block file format
        [microvolts,BlockData] = BlockFileFormat_TOD(fileNames{jj}, EventInfo,jj);
        if ~isempty(BlockData)
            BlockDataAll{1,jj} = BlockData.timestampsMotion;
            BlockDataAll{2,jj} = BlockData.AccelerometerData;
            BlockDataAll{3,jj} = BlockData.GyroscopeData;
            BlockDataAll{4,jj} = BlockData.MagnetometerData;
        else
            BlockDataAll = [];
        end
      end
      
       [filtsig,thetaFiltEEG] = filt_signal(microvolts); % bandpass filter signal
      % Flip signal to invert spikes
        filtsig = -1.*filtsig; 
        thetaFiltEEG = -1.*thetaFiltEEG;
    
    % Common average referencing
      for zz = 1:size(filtsig,1)
        Refsig = median(filtsig(zz,:));
        filtsig(zz,:) = filtsig(zz,:) - Refsig; 
      end
        
    % Separate tetrodes and map DACQ-logger channels
    [tet1_dat,tet2_dat,tet3_dat,tet4_dat] = separate_tets(filtsig); % separate 
    all_tets = {tet1_dat, tet2_dat, tet3_dat, tet4_dat};
    
    [tet1_LFP,tet2_LFP,tet3_LFP,tet4_LFP] = separate_tets(thetaFiltEEG); % separate
    all_LFP{jj,1} = {tet1_LFP, tet2_LFP, tet3_LFP, tet4_LFP};

    % EXTRACT SPIKES % 
    for ii = 1:length(all_tets)
        % Collect spike data
        [spike_mat,spike_count] = extract_spikes(all_tets{ii},ii,timestamps,trialInfo.threshold); 
        % Increment spike count on channel
        num_spikes{ii,jj} = spike_count;   
        % Reshape matrix 
        [final_mat] = reshape_spike_mat(spike_mat);
        % Add spike data to cell array
        new_tets{ii,jj} = final_mat;  
    end
    % Move to next file
    jj = jj + 1;
end
toc

%% Concatenate info from each file
fprintf('#################################################\nConcatenating logger files....\n#################################################\n')
tetrode1 = cell2mat(new_tets(1,:)');
tetrode2 = cell2mat(new_tets(2,:)');
tetrode3 = cell2mat(new_tets(3,:)');
tetrode4 = cell2mat(new_tets(4,:)');
all_LFP = cell2mat(cat(1,all_LFP{:}));

% Sum all spikes in trial
num_spikes = sum(cell2mat(num_spikes),2);

% Motion sensor data: cut to trial time and remove samples when pup touched
if trialInfo.fileFormat == 2 && ~isempty(BlockDataAll)
    BlockDataAll = BlockDataAll';
    motionSensor.msTime = cell2mat(BlockDataAll(:,1));
    motionSensor.accelerometer =  cell2mat(BlockDataAll(:,2));
    motionSensor.gyroscope =  cell2mat(BlockDataAll(:,3));
    motionSensor.magnetometer =  cell2mat(BlockDataAll(:,4));
    
   motionSensor.msTime = motionSensor.msTime - motionSensor.msTime(1);
   
   % Check for duplicate indices
   [~,w] = unique(motionSensor.msTime,'stable');
   duplicate_indices = setdiff( 1:numel(motionSensor.msTime),w)';
   
   structField = fieldnames(motionSensor);
   
    for msInd = 1:size(structField,1)
        motionSensor.(structField{msInd})(duplicate_indices,:) = nan;
        motionSensor.(structField{msInd})(any(isnan(motionSensor.(structField{msInd})), 2), :) = [];
        if msInd == 1      
            trialendInd = find(motionSensor.msTime==(trialInfo.trial_duration*60*1000));            
        end
        motionSensor.(structField{msInd})(trialendInd+1:end,:) = [];   
    end
    
    if ~isempty(EventInfo.TimeEventStart)
        for eventInd = 1:numel(EventInfo.TimeEventStart)
            eventStart = find(motionSensor.msTime==EventInfo.TimeEventStart(eventInd));
            eventEnd = find(motionSensor.msTime==EventInfo.TimeEventEnd(eventInd));          
            for msInd = 2:size(structField,1) % don't remove timestamps
                motionSensor.(structField{msInd})(eventStart:eventEnd,:) = nan;   
            end       
        end
    end
    
else 
    motionSensor = [];
end

%% Cut tetrode matrix to trial time and remove spikes outside this duration
fprintf('Cutting trial to correct length....\n#################################################\n')
[tetrode1,num_spikes(1)] =  cut_trial_to_time(tetrode1,trialInfo,num_spikes(1));
[tetrode2,num_spikes(2)] =  cut_trial_to_time(tetrode2,trialInfo,num_spikes(2));
[tetrode3,num_spikes(3)] =  cut_trial_to_time(tetrode3,trialInfo,num_spikes(3));
[tetrode4,num_spikes(4)] =  cut_trial_to_time(tetrode4,trialInfo,num_spikes(4));


%% Generate Tint .tet header
[genericTetHeader] = TintTetrode_header(trialInfo,DACQinfo);

[tet1FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,1);
[tet2FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,2);
[tet3FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,3);
[tet4FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,4);

%% Append header and create tetrode file
fprintf('Generating tetrode files....\n')
make_tint_file(tetrode1,tet1FileName,1,trialInfo,DACQinfo);
make_tint_file(tetrode2,tet2FileName,2,trialInfo,DACQinfo);
make_tint_file(tetrode3,tet3FileName,3,trialInfo,DACQinfo);
make_tint_file(tetrode4,tet4FileName,4,trialInfo,DACQinfo);

%% Generate EEG file
fprintf('Generating EEG and EGF files....\n#################################################\n')
generateEEG(trialInfo,all_LFP,DACQinfo)

%% Write info to Results struct
Results.trialInfo = trialInfo;
Results.posInfo = posInfo;
Results.tetrode1 = [tetrode1(:,1)./1e6 tetrode1(:,2:end)];
Results.tetrode2 = [tetrode2(:,1)./1e6 tetrode2(:,2:end)];
Results.tetrode3 = [tetrode3(:,1)./1e6 tetrode3(:,2:end)];
Results.tetrode4 = [tetrode4(:,1)./1e6 tetrode4(:,2:end)];
Results.num_spikes = num_spikes;
Results.DACQinfo = DACQinfo;
Results.EventInfo = EventInfo;

% Save Results struct in writeDir
fprintf('Saving Results in MAT-file....\n#################################################\n')
save(fullfile(trialInfo.writeDir,trialInfo.finalTrialName),'Results')
if ~isempty(motionSensor)
    save(fullfile(trialInfo.writeDir,strcat(trialInfo.finalTrialName,'_motionSensor')),'motionSensor')
end

fprintf('END\n') 
end

function [final_mat] = reshape_spike_mat(interp_mat)

interp_mat = permute(interp_mat,[1 3 2]);
final_mat = reshape(interp_mat,[],51); % Reshape to form 2D matrix

end



