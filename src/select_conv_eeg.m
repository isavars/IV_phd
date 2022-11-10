function [bestEEG, peakTheta, peakDelta] = select_conv_eeg(varargin)
%SELECT EEG compares eegs by making power spectrums from a egf files and selects the best one 
%   Takes in data from egf file and outputs bestEEG
%   TO DO: 
%       1. ADD exclude eegs with an snr of under 3 
%       2. speed filter the theta - use speed and position data -
%       each speed point for 5 eeg points - make logical array of indexes
%       (filtering for speeds greater than 2)
tic;
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
 
%load all the eegs and put into a nsamp x 32 matrix called eegs - hc to 32
%but in the future could just find all eegs 

%sleep trial eegs 

sleeptrialname = strcat(data.trials(sleepTrial).trialname,'.eeg');

eegs_S = [];
    for it_eegs = 1:32
        if it_eegs == 1
            eeg_struct  = load_eeg(sleeptrialname);
        else
            eeg_struct = load_eeg(strcat(sleeptrialname,num2str(it_eegs)));
        end
        eegs_S = [eeg_struct.eeg;eegs_S];
    end
    
eegs_S = reshape(eegs_S,[],32); 

%baseline trial eegs

basetrialname = strcat(data.trials(sleepTrial).trialname,'.eeg');

eegs_B = [];
    for it_eegs = 1:32
        if it_eegs == 1
            eeg_struct  = load_eeg(basetrialname);
        else
            eeg_struct = load_eeg(strcat(basetrialname,num2str(it_eegs)));
        end
        eegs_B = [eeg_struct.eeg;eegs_B];
    end
    
eegs_B = reshape(eegs_B,[],32); 

sample_rate = 250;

% get best eeg and delta frequencies from sleep trial: 
    SNR = [];
    for itEEG = 1:32 
        [~, ~, snr] = eeg_powerspec(eegs_S(:,itEEG),sample_rate); 
        SNR(itEEG) = snr;
    end
    [~,bestSNR] = max(SNR);%this finds the index for the eeg with the best SNR 
    peakFreqD = eeg_powerspec(eegs_S(:,bestSNR),sample_rate,'thetaBand',[1.5 4],'hfCutOff',25);

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
    peakFreqT = eeg_powerspec(speedFilteredEEG,sample_rate,'thetaBand',[5 11],'hfCutOff',15);

% outputs 
    bestEEG = eegs_S(:,bestSNR);
    peakTheta = peakFreqT;
    peakDelta = peakFreqD;
toc;
    
end

