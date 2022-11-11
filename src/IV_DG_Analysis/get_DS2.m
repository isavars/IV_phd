function get_DS2 (spatData)
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
    
    sleeptrialname = strcat(trialName{1,6},'.eeg');%check on that 
    
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
    
    [spike_mat,spike_count]=extract_DS2(voltages,sample_rate,spikeThreshold);
    spike_index=0;
    peak_index=10;
    distance_from_peak=2;
    spike_gradient = [];
    while spike_index <= spike_count
        index_of_first_sample_in_cur_batch=spike_index*32+1;
        cur_batch=spike_mat([index_of_first_sample_in_cur_batch:index_of_first_sample_in_cur_batch+32],:);
        pre_spike_col = cur_batch(:,peak_index- distance_from_peak);
        spike_col = cur_batch(:,peak_index);
        post_spike_col = cur_batch(:,peak_index +distance_from_peak);
        pre_grad_col = spike_col - pre_spike_col;
        post_grad_col = post_spike_col - spike_col; 
        spike_gradient = [spike_count, post_grad_col - pre_grad_col]; % spike size made from gradients before and after spike          
    end

    % spike_gradient contains a peak flatness score for 32 channels per
    % spike - should rank flattest at inversion 
        



end