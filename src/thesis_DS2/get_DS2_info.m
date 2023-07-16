function [DS2_info] = get_DS2_info(eeg_data, probe_type, writeDir, sws_table)
%collect the first n DS2s and give mean values to make DS2 info. Then makes and saves DS2 info for
%elePos to be used in classification and DS2 analysis in general. 
% Makes a matrix with DS2 labels per channel for any sleep 
    % trial thats loaded - columns : channel numbers that match spatData 
    % tetrodes, channel depth, DS2 voltages per channel (amplitude), DS2 label
    % per channel ("above", "in", "below"), distance from DS2 inversion in
    % absolute depth (only possible on shanks with inversion) 
% INPUT - eeg_data - file path to egf file from sleep trial, probe_type - 1
% single opp, 2 - multi opp, 3 - single same, 
%TO DO 
% 1) add eeg_channel identities to DS2 so eegs can be matched to tetrodes
% later in analysis. 

    %get trial name for saving files (probably there's a better way to do
    %this) 
    lastBackslashIdx = find(eeg_data == '\', 1, 'last');
    sleep_eeg = eeg_data(lastBackslashIdx+1:end);
    sleep_trial = extractBefore(sleep_eeg, '.egf');
    %load SWS data from table to get parameters 
    load(sws_table, 'sleepData');  
    curr_trial = find(sleepData.trialname == sleep_trial);
    dataset = char(sleepData.dataset(curr_trial));
    
    %call get_DS2 to get all the necessary bits 
 
    [max_amplitude, mean_amplitude, fig1, fig2, ds2_rate, peak_to_trough_amplitude, peak_to_trough_slope, spike_mat, new_eeg_chans] = get_DS2_from_egf( eeg_data, probe_type,sws_table);

    %save spike_mat for future data processing and 1st Ds2 
    %make sure you've made the subfolders 
    if ~isempty(spike_mat)
        save ([writeDir '\spike_mats\' char(dataset) '_spike_mat.mat'], 'spike_mat')
        savefig(fig1, [writeDir '\figures\' char(dataset) '_DS2_spikes_figure.fig'])
        try
            savefig(fig2, [writeDir '\figures\' char(dataset) '_DS2_spike_vs_depth_figure.fig'])
        catch ME 

        end
    else
    end

    %produce DS2 labels (in a table) per channel per rat per day-  makes 3 kinds of table depending on
    %probe type. 

    % Define the variable names and types
    varNames = {'channel', 'channel_depth', 'DS2_labels', 'DS2_spike_times', 'DS2_max_amplitude', 'DS2_mean_amplitude', 'DS2_rate','DS2_peak_to_trough_amplitude','DS2_slope'}; %'DS2_distance'}; %temporary place for DS2 rate 
    varTypes = {'double', 'double', 'string', 'cell', 'double', 'double', 'double','double', 'double'};

    % Create an empty table with the specified shape
    DS2_info = table('Size', [32, length(varTypes)], 'VariableNames', varNames, 'VariableTypes', varTypes);
    %add channel numbers to the table (same for all silicone probes- changes below for tetrodes) 
    DS2_info.channel = (1:32)'; 
    % add ds2_rate along every channel - not needed but will transfer over
    % in a different way later 
    DS2_info.DS2_rate = repmat(ds2_rate, 32, 1);
    spike_times = cell(size(spike_mat,3),1);
    %add DS2 spiketimes - should be the same for every probetype 
    for ii = 1: size(spike_mat,3) %iterate over spikes and get all the times
        spike_times{ii} = spike_mat(1,1,ii); %1 on first dimension is because you only need one of the channels, 1 on second dimension is where the spike time is kept on the array.          
    end
    for jj = 1:32
        DS2_info.DS2_spike_times{jj} = spike_times;
    end 
    %get the maximum max amplitude so ds2 labels can be made as a
    %proporiton of the size of the largest ds2 spike from that rat 
    max_max_amplitude = nanmax(abs(max_amplitude));
    
    % insert data into the table 
    if probe_type == 1 
        load('single_shank_probe_map.mat','single_shank_probe_map')
        for it_ch = 1:32 
            DS2_info.channel_depth(it_ch) = single_shank_probe_map(it_ch,3);%check that mapping makes sense - does this need to be converted? 
            %make features from selected spike - this big should be a
            %function that gets called once per probe 
            DS2_info.DS2_max_amplitude(it_ch) = max_amplitude(it_ch);
            DS2_info.DS2_mean_amplitude(it_ch) = mean_amplitude(it_ch);
            DS2_info.DS2_peak_to_trough_amplitude(it_ch) = peak_to_trough_amplitude(it_ch);
            DS2_info.DS2_slope(it_ch) = peak_to_trough_slope(it_ch);
            if max_amplitude(it_ch) < - abs(max_max_amplitude*0.5) 
                DS2_info.DS2_labels(it_ch) = "above";
            elseif max_amplitude(it_ch) >= - abs(max_max_amplitude*0.5) && max_amplitude(it_ch) <= abs(max_max_amplitude*0.5)
                DS2_info.DS2_labels(it_ch) = "in";
            elseif max_amplitude(it_ch) > abs(max_max_amplitude*0.5)
                DS2_info.DS2_labels(it_ch) = "bellow";
            end
        end
    elseif probe_type == 2
        load('multi_shank_probe_map.mat','multi_shank_probe_map')
        for it_ch = 1:32              
            DS2_info.channel_depth(it_ch) = multi_shank_probe_map(it_ch,3); 
            DS2_info.DS2_max_amplitude(it_ch) = max_amplitude(it_ch); % fill amplitudes per channel just by indexing the voltages 
            DS2_info.DS2_mean_amplitude(it_ch) = mean_amplitude(it_ch);
            DS2_info.DS2_peak_to_trough_amplitude(it_ch) = peak_to_trough_amplitude(it_ch);
            DS2_info.DS2_slope(it_ch) = peak_to_trough_slope(it_ch);
            if max_amplitude(it_ch) < - abs(max_max_amplitude*0.5) 
                DS2_info.DS2_labels(it_ch) = "above";
            elseif max_amplitude(it_ch) >= - abs(max_max_amplitude*0.5) && max_amplitude(it_ch) <= abs(max_max_amplitude*0.5)
                DS2_info.DS2_labels(it_ch) = "in";
            elseif max_amplitude(it_ch) > abs(max_max_amplitude*0.5)
                DS2_info.DS2_labels(it_ch) = "bellow";
            end
            % DS2_info.DS2_distance = can use spread value to get max and
            % ming then chop into thirds to create the three possible
            % distances from inversion 
            %will need to label the inflection point 0 um and then only label
            % chanels on shanks with inflection points with distance from
            % that 0 in um using 'channel_depth'
        end
    elseif probe_type ==3 
        load('single_shank_same_probe_map_final.mat','single_shank_probe_map')
        for it_ch = 1:32
            DS2_info.channel_depth(it_ch) = single_shank_probe_map(it_ch,3);
            %make amplitude and label features from selected spike 
            DS2_info.DS2_max_amplitude(it_ch) = max_amplitude(it_ch);
            DS2_info.DS2_mean_amplitude(it_ch) = mean_amplitude(it_ch);
            DS2_info.DS2_peak_to_trough_amplitude(it_ch) = peak_to_trough_amplitude(it_ch);
            DS2_info.DS2_slope(it_ch) = peak_to_trough_slope(it_ch);
            if max_amplitude(it_ch) < - abs(max_max_amplitude*0.5) 
                DS2_info.DS2_labels(it_ch) = "above";
            elseif max_amplitude(it_ch) >= - abs(max_max_amplitude*0.5) && max_amplitude(it_ch) <= abs(max_max_amplitude*0.5)
                DS2_info.DS2_labels(it_ch) = "in";
            elseif max_amplitude(it_ch) > abs(max_max_amplitude*0.5)
                DS2_info.DS2_labels(it_ch) = "bellow";
            end
        end
    elseif probe_type ==4 
        %load channel length 
        num_eegs = size(sleepData.eeg_channels{curr_trial},1);
%         eeg_chans = sleepData.eeg_channels{curr_trial}; %using
%         new_eeg_chans intead. 

        for it_ch = 1:32
            %no channel depth info is availabe for tetrodes 
            DS2_info.channel_depth(it_ch) = nan;
            %change channel labels to correspond to each signal
            if it_ch  <= length(new_eeg_chans)
                DS2_info.channel(it_ch) = new_eeg_chans(it_ch);
            else 
                DS2_info.channel(it_ch) = nan;
            end
            %make amplitude and label features from selected spike 
            DS2_info.DS2_max_amplitude(it_ch) = max_amplitude(it_ch);
            DS2_info.DS2_mean_amplitude(it_ch) = mean_amplitude(it_ch);
            DS2_info.DS2_peak_to_trough_amplitude(it_ch) = peak_to_trough_amplitude(it_ch);
            DS2_info.DS2_slope(it_ch) = peak_to_trough_slope(it_ch);
            if max_amplitude(it_ch) < - abs(max_max_amplitude*0.5) 
                DS2_info.DS2_labels(it_ch) = "above";
            elseif max_amplitude(it_ch) >= - abs(max_max_amplitude*0.5) && max_amplitude(it_ch) <= abs(max_max_amplitude*0.5)
                DS2_info.DS2_labels(it_ch) = "in";
            elseif max_amplitude(it_ch) > abs(max_max_amplitude*0.5)
                DS2_info.DS2_labels(it_ch) = "bellow";
            end
        end
    end

    
    save ([writeDir '\' char(dataset) '_DS2_info.mat'], 'DS2_info')
    

end 