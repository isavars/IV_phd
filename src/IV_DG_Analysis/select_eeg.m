function [bestEEG, peakTheta, peakDelta] = select_eeg(varargin)
%SELECT EEG compares eegs by making power spectrums from a scan file and selects the best one 
%   Takes in data from scan file and outputs bestEEG
%   TO DO: 
%       1. ADD exclude eegs with an snr of under 3 
%       2. speed filter the theta - use speed and position data -
%       each speed point for 5 eeg points - make logical array of indexes
%       (filtering for speeds greater than 2)

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
    SNR = [];
    for itEEG = 1:length(data.trials(sleepTrial).eeg)
        [~, ~, snr] = eeg_powerspec(data.trials(sleepTrial).eeg(itEEG).eeg,data.trials(sleepTrial).eeg(itEEG).sample_rate); 
        SNR(itEEG) = snr;
    end
    [~,bestSNR] = max(SNR);%this finds the index for the eeg with the best SNR 
    peakFreqD = eeg_powerspec(data.trials(sleepTrial).eeg(bestSNR).eeg,data.trials(sleepTrial).eeg(bestSNR).sample_rate,'thetaBand',[1.5 4],'hfCutOff',25);

% get peak theta from baseline trial:
    speedFilteredEEG = [];
    baselineEEG = data.trials(baselineTrial).eeg(bestSNR).eeg;
    filteredSpeed = data.trials(baselineTrial).speed >= 2;% this needs to account for the 5 to 1 eeg to speed 
    filteredSpeed = reshape(repmat(filteredSpeed,1,5).',1,[]);
    for itSp = 1: length(filteredSpeed)
         if filteredSpeed(itSp) == 1 
            speedFilteredEEG = [speedFilteredEEG; baselineEEG(itSp)];
         end
    end
    peakFreqT = eeg_powerspec(speedFilteredEEG,data.trials(baselineTrial).eeg(bestSNR).sample_rate,'thetaBand',[5 11],'hfCutOff',15);

% outputs 
    bestEEG = data.trials(sleepTrial).eeg(bestSNR).eeg;
    peakTheta = peakFreqT;
    peakDelta = peakFreqD;
    
end

