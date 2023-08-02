function spatData = make_wf_and_sleep_features(read_dir, spatial_data, sleep_data)
%this function loops through the all waveform matrices per
%experiment and makes features from the waveform data that can be loaded into
%spatData variables so they are the length of spatData. 
%AND adds sleep information using the sleepData table - this should add 
%rem info - duration, and spike times, mean rate and nspks for rem part of 
% trial. It should also modify the spkTs, meanRate, peakrate, nSPks and
% burst index for the sleep trial column and the trial duration for the
% sleep trial also 

% ISSUE #1 - datasets in spatData need to be in acensding order so unique
% function works as an index.
%#2 needs to work with tetrode data so 


% Inputs - read_dir (path to folder with experiment based waveform
% matrices, spatial_data (spatData), sleep_data (sleepData)
% Outputs - waveforms (max channel, mean wf 97x1), wf_means (mean of all
% the spikes in the cluster (97x4), trough-to-peak latency (uses waveform 
% from waveforms - lowest value (after peak) - highest value should be trough peak).

    % load spatData 
    load (spatial_data, 'spatData')

    %get all the waveform features to be added to spatData
    [waveforms, wf_means, max_wf_channel, TP_latency] = make_wf_features(read_dir, spatData); 


    %load sleepData
    load (sleep_data, 'sleepData')

    %get all the sleep features so they can be added to spatData 
    
    [swsTrialDuration, swsSpikeTimes, remTrialDuration, remSpikeTimes, swsMeanRate, swsNSpks, swsBurstIndex , remMeanRate, remNSpks, remBurstIndex] = make_sleep_features(spatData, sleepData);

    %add new features to spatData and save 

    %waveform features 
    spatData.waveforms = waveforms; 
    spatData.wf_means = wf_means;
    spatData.max_wf_channel = max_wf_channel;
    spatData.TP_latency = TP_latency;
    
    %sleep modifications
    % get sleep trial index and add to the correct part of the table 
    for ii= 1:height(spatData)
        sleep_trial = strcmp(string(spatData.env(ii,:)),'sleep');
        sleep_idx = find(sleep_trial,1);
%         spatData.trialDur(ii,sleep_idx) = swsTrialDuration(ii); 
        spatData.SpkTs(ii,sleep_idx) = swsSpikeTimes(ii);
        spatData.meanRate(ii,sleep_idx) = swsMeanRate(ii); % mean rate 
        spatData.nSpks(ii,sleep_idx) = swsNSpks(ii); % number of spikes    
        spatData.burstIndex(ii,sleep_idx) = swsBurstIndex(ii);
    end
    
    %rem sleep info and sws trial duration
    spatData.swsTrialDur = swsTrialDuration;

    spatData.remTrialDur = remTrialDuration; 
    spatData.remSpkTs = remSpikeTimes; %make sure these are a cell array 
    spatData.remMeanRate = remMeanRate; % mean rate 
    spatData.remNSpk = remNSpks; % number of spikes    
    spatData.remBurstIndex = remBurstIndex;


%     save(spatial_data,"spatData");

end

function [waveforms, wf_means, max_wf_channel, TP_latency] = make_wf_features(read_dir, spatData)
    % make index for single experiments on spatData 
    unique_dataset = unique(spatData.dataset); % rearagne to match spatData order 
    spatData_idx = [];
    for itS = 1: height(spatData)
        spatData_idx(itS) = find(strcmp(unique_dataset, spatData.dataset(itS)));
    end
    spatData_idx = spatData_idx';

    %sampling rate to be used for trough-to-peak computation 
    sample_rate = 48; %in kHz so the converion to ms is more straightforward 

    %make the output features
    waveforms = [];
    wf_means  = [];
    max_wf_channel = [];
    TP_latency= [];

    % loop throught all datasets on spatData and get waveforms and wf feaures. 
    for itD = 1:(length(unique_dataset))
        %load correct dataset
        filename = fullfile(read_dir, [unique_dataset{itD} '_wfs.mat']);
        %make indexes for dataset within spatData -pretty sure it should be
        %here in the loop 
        data_idx = find(spatData_idx == itD, 1,'first');
        data_idx_2 = find(spatData_idx == itD, 1,'last');
        %make temporary waveform features 
        waves = cell(length(data_idx:data_idx_2), (size(spatData.trialNo,2)) ); 
        max_wf = cell(length(data_idx:data_idx_2), (size(spatData.trialNo,2)) );
        wf_chan = zeros(length(data_idx:data_idx_2), (size(spatData.trialNo,2)) );
        temp_TPL = zeros(length(data_idx:data_idx_2), (size(spatData.trialNo,2)) );

        spatData_curr_trial = spatData(data_idx:data_idx_2, :);
               
        %loop through trial
        for trial_it = 1: sum(~isnan(spatData.trialNo(data_idx,:))) %adjusted for varying trial number lengths
            %loop through cell
            for cell_it = 1: length(data_idx:data_idx_2)
                if isfile(filename)
                    load(filename,'exp_wfs');
                    %make waveform means
                    waves{cell_it,trial_it} = mean(exp_wfs{cell_it,trial_it},3);
                else %this is unecessary since the values are already in the table but adding so it runs as is              
                    waves{cell_it,trial_it} = spatData_curr_trial.wf_means{cell_it,trial_it};
                end
                %find the mean waveform for the channel with maximum
                %amplitude 
                 if ~isempty(waves{cell_it,trial_it})
                    temp_wf = waves{cell_it,trial_it};
                    wfMins = nanmin(temp_wf, [], 1 );
                    wfMaxs = nanmax(temp_wf, [], 1 );
                    wfAmps = wfMaxs - wfMins;
                    [~,maxAmpCh] = nanmax( wfAmps );
                    wf_chan(cell_it, trial_it) = maxAmpCh;
                    best_wf = temp_wf( :, maxAmpCh );
                    max_wf{cell_it, trial_it} = best_wf;
                    %calculate trough-to-peak latency per max waveform
                    peak_amp = max(best_wf);
                    peak_amp_idx = find(peak_amp == best_wf);
                    trough_amp = min(best_wf(peak_amp_idx:end)); %min needs to be after the peak
                    peak_time = (find(peak_amp == best_wf,1))./sample_rate;
                    trough_time = ((find(trough_amp == best_wf,1)))./sample_rate;               
                    temp_TPL(cell_it, trial_it) = trough_time - peak_time;
                else 
                    max_wf{cell_it,trial_it} = nan;
                    wf_chan(cell_it,trial_it) = nan;
                    temp_TPL(cell_it,trial_it) = nan;
                 end 
            end        
        end 
        wf_means = [wf_means; waves];
        waveforms = [waveforms; max_wf];
        max_wf_channel = [max_wf_channel; wf_chan];
        TP_latency = [TP_latency; temp_TPL];
        clear exp_wfs  
    end
end
function [swsTrialDuration, sws_SpikeTimes, remTrialDuration, rem_SpikeTimes, swsMeanRate, swsNSpks, swsBurstIndex, remMeanRate, remNSpks, remBurstIndex] = make_sleep_features(spatData, sleepData)
    %for every row of spatData - if a dataset is in sleepData use the
    %sws_epoch and rem_epochs from that dataset to make new spike times and
    %new trial duraitons for each type of sleep   
    
    sws_SpikeTimes = cell(height(spatData),1);
    rem_SpikeTimes = cell(height(spatData),1);
    swsTrialDuration =[];
    remTrialDuration =[];
    %features made from spike times and trial durations 
    swsMeanRate = [];
    swsNSpks = [];
    swsBurstIndex = [];

    remMeanRate = [];
    remNSpks = [];
    remBurstIndex = [];


    for itS = 1: height(spatData)
       sleepData_idx = find(strcmp(spatData.dataset(itS), sleepData.dataset),1);
       sws_epochInds = sleepData.SWS_epochs{sleepData_idx};
       rem_epochInds = sleepData.REM_epochs{sleepData_idx};

       %find index for sleep trial 
       sleep_trial = strcmp(string(spatData.env(itS,:)),'sleep');
       sleep_idx = find(sleep_trial,2);
       if size(sleep_idx,2) > 1 
           sleep_idx = sleep_idx(2); %dealing with trial with more than one sleep
       end
       spikeTimes = spatData.SpkTs{itS,sleep_idx}; %spike times per cell for sleep trial (to be cut) 
    
       %make durations per dataset 
       swsTrialDuration = [swsTrialDuration; sleepData.duration(sleepData_idx)];  
       remTrialDuration = [remTrialDuration; sleepData.REM_duration(sleepData_idx)];

       %make spiketimes per dataset 
       swsSpikeTimes =[];
       for itEps=1:size(sws_epochInds, 2)
           indEps = spikeTimes >= sws_epochInds(1,itEps) & spikeTimes <= sws_epochInds(2,itEps);
           swsSpikeTimes = [swsSpikeTimes; spikeTimes(indEps,:)];   
       end
       remSpikeTimes = [];
       for itEps=1:size(rem_epochInds, 2)
           indEps = spikeTimes >= rem_epochInds(1,itEps) & spikeTimes <= rem_epochInds(2,itEps);
           remSpikeTimes = [remSpikeTimes; spikeTimes(indEps,:)];   
       end   

       swsMeanRate = [swsMeanRate; length(swsSpikeTimes) / swsTrialDuration(1)]; % mean rate 
       swsNSpks = [swsNSpks; length(swsSpikeTimes )]; % number of spikes    
       swsBurstIndex = [swsBurstIndex; (sum(diff(swsSpikeTimes) <= 0.009))/(length(diff(swsSpikeTimes)))];
    
       remMeanRate = [remMeanRate; length(remSpikeTimes) / remTrialDuration(1)]; % mean rate 
       remNSpks = [remNSpks; length(remSpikeTimes )]; % number of spikes    
       remBurstIndex = [remBurstIndex; (sum(diff(remSpikeTimes) <= 0.009))/(length(diff(remSpikeTimes)))]; 

       sws_SpikeTimes{itS} = swsSpikeTimes;
       rem_SpikeTimes{itS} = remSpikeTimes;


    end   
   
end 