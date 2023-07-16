function [sleepData] = getSleepData() 
%this function gets sleep data from trials loded into a scan and saves them
%into a table to be used by either getSpatData , or get_DS2 (doing it just
%for DS2 now) 
%INPUT - path to write directory OUTPUT -table with sws and rem epochs save it in a
%desired location when you are done 
%TO DO 1) add pathname to write directory and save as 'sleepData' 2)can
%make .mats of eegs only make for all the relevant days in between trials -
%make sure this function runs on one of those 

SD=gss;%grab all scan datasets 

 %make table sleepData
    varNames = {'ratID', 'dataset', 'trialname','trialpath' ,'eeg_channels','duration', 'SWS_epochs','REM_duration', 'REM_epochs'  };
    varTypes = {'string','string','string', 'string', 'cell', 'double', 'cell', 'double', 'cell'}; %sws_inds will have two columns of differnt legths per trial
    
    sleepData = table('Size', [length(SD.selData), length(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);

%make loop to collect data from each file loaded on scan 
    for ii = 1:length(SD.selData)
           data = evalin('base', SD.selData{ii}); %load current dataset        
           [bestEEG, peakTheta, peakDelta,sleepTrial] = select_eeg(SD,ii); 
           stateInds =  getBrainStateHardThr (data.trials(sleepTrial).speed, bestEEG, 'thetaFr', peakTheta, 'deltaFr', peakDelta);
           %get SWS indices to cut the trials down to only sleep 
           SWS = stateInds.sws;
           indStart = diff([0,SWS,0]) == 1;  indStart = find(indStart == 1); indStart = indStart*0.8 - 0.4;
           indStop = diff([0,SWS,0]) == -1; indStop = find(indStop == 1); indStop = indStop*0.8 + 0.4; %what do 0.8 and 0.4 mean here?
           SWS_epochs = [indStart; indStop]; 
           %filter out any sleep epochs that are shorter than 5s to reduce
           %noise
           for it_eps = 1:length(SWS_epochs)
                if diff(SWS_epochs(:,it_eps)) < 5
                    SWS_epochs(:,it_eps) = nan;
                else
                end
           end
           SWS_epochs = reshape(SWS_epochs(~isnan(SWS_epochs)),2, []); %reshape into SWS_epochs 
           newTrialDuration = sum(diff(SWS_epochs));
           %get REM indices incase you want to use them for CA3 filtering 
           REM = stateInds.rem;
           indStart = diff([0,REM,0]) == 1;  indStart = find(indStart == 1); indStart = indStart*0.8 - 0.4;
           indStop = diff([0,REM,0]) == -1; indStop = find(indStop == 1); indStop = indStop*0.8 + 0.4; %what do 0.8 and 0.4 mean here?
           REM_epochs = [indStart; indStop]; 
           remTrialDuration = sum(diff(REM_epochs));

           %get channel numbers for eggs - trials is 1 here because same
           %eegs were recorded in every trial. 
           eeg_channels = cell(length(data.trials(1).eeg),1);
           for jj = 1: length(data.trials(1).eeg)
                eeg_channels{jj} = data.trials(1).eeg(jj).channel;
           end 

           %add values to the table 
           sleepData.ratID(ii) = cellstr(data.user.ID);
           sleepData.dataset(ii) = cellstr(data.load_selection.dataName);
           sleepData.trialname(ii) = data.trials(sleepTrial).trialname; 
           sleepData.trialpath(ii) =  cellstr(data.load_selection.pathName); 
           sleepData.eeg_channels{ii} = eeg_channels; %I think this shoudl work 
           sleepData.duration(ii) = newTrialDuration;
           sleepData.SWS_epochs{ii} = SWS_epochs;
           %add REM
           sleepData.REM_duration(ii) = remTrialDuration;
           sleepData.REM_epochs{ii} = REM_epochs;
           
    end 

end

function [bestEEG, peakTheta, peakDelta, sleepTrial] = select_eeg(SD,exp_idx)
%SELECT EEG compares eegs by making power spectrums from a scan file and selects the best one 
%   Takes in data from scan file and outputs bestEEG
%   TO DO: 
%       1. ADD exclude eegs with an snr of under 3 
%       2. speed filter the theta - use speed and position data -
%       each speed point for 5 eeg points - make logical array of indexes
%       (filtering for speeds greater than 2) %01/06/23 Iv - is this fixed?       

    %get current selecter dataset in scan 
    data = evalin('base', SD.selData{exp_idx});
    for  itTR = 1:length(data.trials) % make indexes for baseline and sleep 
        if strcmp(data.trials(itTR).user.environment, 'fam')
            baselineTrial = itTR;
        elseif strcmp(data.trials(itTR).user.environment, 'sleep')
            sleepTrial = itTR;
        end
    end
    
% get best eeg and delta frequencies from sleep trial: 

    %get eegs for sleep trial and baseline trial (eegs need to be in the file the scan .mat was made from) 
    sleeptrialname = [data.load_selection.pathName data.trials(sleepTrial).trialname '.eeg'];
    baselinetrialname = [data.load_selection.pathName data.trials(baselineTrial).trialname '.eeg'];
    
    %find number of files containing the sleeptrialname 
    fileList = dir(fullfile(data.load_selection.pathName, '*.eeg*'));
    eeg_count = 0;
    eeg_ext = [];
    % Iterate through each file in the folder
    for i = 1:numel(fileList)
        % Extract the file name
        [~, fileName, fileExt] = fileparts(fileList(i).name);        
        % Check if the file name contains the search string
        if contains(string([fileName '.eeg']), [data.trials(sleepTrial).trialname '.eeg'])
            eeg_count = eeg_count + 1;
            file_ext = str2double(extractAfter( fileExt, '.eeg'));
            if isempty(file_ext)
                eeg_ext = [eeg_ext; 1];
            else 
                eeg_ext = [eeg_ext; file_ext];
            end
        end
    end

    eegs_S = [];
    eegs_B = [];
    for it_eegs = 1:eeg_count
        if it_eegs == 1
            eeg_struct_s  = load_eeg(sleeptrialname);
            eeg_struct_b  = load_eeg(baselinetrialname);
        else
            eeg_struct_s = load_eeg(strcat(sleeptrialname,num2str(eeg_ext(it_eegs))));
            eeg_struct_b = load_eeg(strcat(baselinetrialname,num2str(eeg_ext(it_eegs))));
        end
        eegs_S = [eeg_struct_s.eeg;eegs_S];
        eegs_B = [eeg_struct_b.eeg;eegs_B];
    end
    eegs_S = reshape(eegs_S,[],eeg_count);
    eegs_B = reshape(eegs_B,[],eeg_count);

    SNR = [];    
    for itEEG = 1:eeg_count 
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