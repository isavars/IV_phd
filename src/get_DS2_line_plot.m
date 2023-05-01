function get_DS2_line_plot(spatData, spike_idx)
% get_DS2 takes eegs from all chanels and filters them to detect DS2 spikes
% it then compares DS2 spike amplitude across channels to detect an
% inversion point (0 amp) and labels all channels at/around the point and
% above GCL and all below HL 

    %load all the eegs and put into a nsamp x 32 matrix called eegs - hc to 32
    %but in the future could just find all eegs - 
    
    %sleep trial eegs 
    
    %find the sleep trial using info in spatData - add trial labels to the
    %table
    trialName = spatData.trialName;
    
%     sleeptrialname = strcat(trialName{1,6},'.eeg');%check on that 

    sleeptrialname = '220505j_sleepHP.eeg';
    
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
    
    %convert to mV 
    voltages = eegs_S.*(1.5 ./ 128);
    voltages = voltages'; 
    
    %run extract spikes TW but put in Senzai and Buzsaki parameters for DS2
    %spikes 

    sample_rate = 250;
    spikeThreshold = 1.14;%based on Senzai and Buzsaki Threshold 
    

    % extract DS2 from eegs 
    [spike_mat,spike_count]=extract_DS2_TV(voltages,sample_rate,spikeThreshold);
    %[~, max_spread_idx] = max(spike_mat(1, 2, :),[],'all', 'linear');
    [~, spike_idx_sorted_by_spread] = sort(spike_mat(1, 2, :), 'descend');

    % Hyperparameter here to only look at top 100 events by spread
    top_k_spikes = 100;
    top_k_mat = nan(top_k_spikes, 2, 32);
    for k = 1:top_k_spikes
%         spike_idx = spike_idx_sorted_by_spread(k); %commented out so it runs on selected spike 
        probe_data = load('multi_shank_probe_map.mat');%load('single_shank_probe_map.mat');
        probe_data = struct2array(probe_data);
        probe_idx = probe_data(:, 1) + 1;
        probe_depth = probe_data(:, 3);

        for i = 1:length(probe_idx)
            p_idx = probe_idx(i);
            top_k_mat(k, 1, i) = spike_mat(p_idx,20,spike_idx);
            top_k_mat(k, 2, i) = probe_depth(i);
        end
    end

    %makes a line plot of a single event with voltage on the x and depth on the y  
    
%     selected_spike = top_k_mat(spike_idx,:,:);
%     selected_spike = squeeze(selected_spike);
%     
%     figure()
%     plot(selected_spike(1,:),selected_spike(2,:));

    % The next step would be to take what we've done above and repeat it
    % for K candidate events. We would then average the per depth voltage
    % AT putative DS2 spike.

    figure()
    avg_mat = squeeze(mean(top_k_mat, 1));
    whos
    plot(avg_mat(1, :), avg_mat(2, :));


    % we want to get the max, per spike window and its t_max
    % we want to get the value for t_max+1





end