function [waveforms, wf_means, TP_latency] = make_wf_features(read_dir, data)
%this function loops through the all waveform matrices per
%experiment and makes features from the waveform data that can be loaded into
%spatData variables so they are the length of spatData. 
% ISSUE #1 - datasets in spatData need to be in acensding order so unique
% function works as an index.
% ISSUE #2 - shouldn't be 5 trials wide - need to get sleep trial data too
% or at least allow space for it. 

% Inputs - read_dir (path to folder with experiment based waveform
% matrices, data (spatData).
% Outputs - waveforms (max channel, mean wf 97x1), wf_means (mean of all
% the spikes in the cluster (97x4), trough-to-peak latency (uses waveform 
% from waveforms - lowest value (after peak) - highest value should be trough peak).

    % load spatData 
    load (data, 'spatData')

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

    % loop throught all datasets on spatData and get mean waveforms 
    for itD = 1:(length(unique_dataset))
        %load correct dataset
        load(fullfile(read_dir, [unique_dataset{itD} '_wfs.mat']),'exp_wfs')
            %disp(exp_wfs)
        %make indexes for dataset within spatData
        data_idx = find(spatData_idx == itD, 1,'first');
        data_idx_2 = find(spatData_idx == itD, 1,'last');
        %make temporary waveform features 
        waves = cell(length(data_idx:data_idx_2), (size(spatData.trialNo,2)) -1); %will need to edit these to it grabs the sleep trial 
        max_wf = cell(length(data_idx:data_idx_2), (size(spatData.trialNo,2)) -1);
        wf_chan = zeros(length(data_idx:data_idx_2), (size(spatData.trialNo,2)) -1);
        temp_TPL = zeros(length(data_idx:data_idx_2), (size(spatData.trialNo,2)) -1);
        %loop through trial
        for trial_it = 1: (size(spatData.trialNo,2)) -1
            %loop through cell
            for cell_it = 1: length(data_idx:data_idx_2)
                %make waveform means
                waves{cell_it,trial_it} = mean(exp_wfs{cell_it,trial_it},3);
                %find the mean waveform for the channel wiht maximum
                %amplitude 
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
            end        
        end 
        wf_means = [wf_means; waves];
        waveforms = [waveforms; max_wf];
        max_wf_channel = [max_wf_channel; wf_chan];
        TP_latency = [TP_latency; temp_TPL];
        clear exp_wfs
    end 
    
    %add new features to spatData and save 
    spatData.waveforms = waveforms; 
    spatData.wf_means = wf_means;
    spatData.max_wf_channel = max_wf_channel;
    spatData.TP_latency = TP_latency;

    save(data,"spatData");

end