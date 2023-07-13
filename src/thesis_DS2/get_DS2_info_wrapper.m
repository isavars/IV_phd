function [DS2_info] = get_DS2_info_wrapper(num_spikes, eeg_data, probe_type, writeDir, sws_table)
%collect the first n DS2s and give mean values to make DS2 info. Then makes and saves DS2 info for
%elePos to be used in classification and DS2 analysis in general. 
% Makes a matrix with DS2 labels per channel for any sleep 
    % trial thats loaded - columns : channel numbers that match spatData 
    % tetrodes, channel depth, DS2 voltages per channel (amplitude), DS2 label
    % per channel ("above", "in", "below"), distance from DS2 inversion in
    % absolute depth (only possible on shanks with inversion) 
%TO DO 1)  reducde the number of input arguments - numspikes needs to be 
% inside of get_DS2_from_egf. 2) combine with run_ds2 wrappers 3) doesnt
% need to make one DS2_info table per probe unless you want channel depth
% and channel depth can be coded into the pipeline somewhere else since
% probe type is a per rat feature regardless. 

    %get trial name for saving files (probably there's a better way to do
    %this) 
    lastBackslashIdx = find(eeg_data == '\', 1, 'last');
    sleep_eeg = eeg_data(lastBackslashIdx+1:end);
    sleep_trial = extractBefore(sleep_eeg, '.eeg');
    %load SWS data from table to get parameters 
    load(sws_table, 'sleepData');  
    curr_trial = find(sleepData.trialname == sleep_trial);
    dataset = char(sleepData.dataset(curr_trial));
    
    %call get_DS2 to get all the necessary bits - spike_mat is made here
    %once per spike_num this is redundant. spike_num needs to itterate over
    %a spike_mat made once - all the info is in spike_mat right? 
    amplitudes = [];
    mean_amplitudes = [];
    for spike_num = 1:num_spikes    
        [amplitude, mean_amplitude, fig1, fig2, ds2_rate, peak_to_trough_amplitude, peak_to_trough_slope, spike_mat] = get_DS2_from_egf( eeg_data,spike_num, probe_type,sws_table,curr_trial);
        amplitudes = [amplitudes, amplitude];
        mean_amplitudes = [mean_amplitudes, mean_amplitude];
        cur_spike = num2str(spike_num);
        %savefig ([writeDir '\' rat_ID '_' sleep_trial '_DS2_spike_figure_' cur_spike '.mat'],'fig1','fig2')
    end 
    amplitude = mean(amplitudes,2);
    mean_amplitude = mean(mean_amplitudes,2);

    %save spike_mat for future data processing (commented out until I'm
    %making final version of elePos)

    %save ([writeDir '\' char(dataset) '_spike_mat.mat'], 'spike_mat')

    %produce DS2 labels (in a table) per channel per rat per day-  makes 3 kinds of table depending on
    %probe type. 

    % Define the variable names and types
    varNames = {'channel', 'channel_depth', 'DS2_labels', 'DS2_max_amplitude', 'DS2_mean_amplitude', 'DS2_rate','DS2_peak_to_trough_amplitude','DS2_slope'}; %'DS2_distance'}; %temporary place for DS2 rate 
    varTypes = {'double', 'double', 'string', 'double', 'double', 'double','double', 'double'};

    % Create an empty table with the specified shape
    DS2_info = table('Size', [32, 6], 'VariableNames', varNames, 'VariableTypes', varTypes);
    %add channel numbers to the table (same for all probes) 
    DS2_info.channel = (1:32)'; 
    % add ds2_rate along every channel - not needed but will transfer over
    % in a different way later 
    DS2_info.DS2_rate = repmat(ds2_rate, 32, 1);
    
    % insert data into the table 
    if probe_type == 1 
        load('single_shank_probe_map.mat','single_shank_probe_map')
        for it_ch = 1:32 
            DS2_info.channel_depth(it_ch) = single_shank_probe_map(it_ch,3);%check that mapping makes sense - does this need to be converted? 
            %make features from selected spike 
            DS2_info.DS2_max_amplitude(it_ch) = amplitude(it_ch);
            DS2_info.DS2_mean_amplitude(it_ch) = mean_amplitude(it_ch);
            DS2_info.DS2_peak_to_trough_amplitude(it_ch) = peak_to_trough_amplitude(it_ch);
            DS2_info.DS2_slope(it_ch) = peak_to_trough_slope(it_ch);
            if amplitude(it_ch) < -0.5 
                DS2_info.DS2_labels(it_ch) = "above";
            elseif amplitude(it_ch) >= -0.5 && amplitude(it_ch) <= 0.5
                DS2_info.DS2_labels(it_ch) = "in";
            elseif amplitude(it_ch) > 0.5
                DS2_info.DS2_labels(it_ch) = "bellow";
            end
        end
    elseif probe_type == 2
        load('multi_shank_probe_map.mat','multi_shank_probe_map')
        for it_ch = 1:32              
            DS2_info.channel_depth(it_ch) = multi_shank_probe_map(it_ch,3); 
            DS2_info.DS2_max_amplitude(it_ch) = amplitude(it_ch); % fill amplitudes per channel just by indexing the voltages 
            DS2_info.DS2_mean_amplitude(it_ch) = mean_amplitude(it_ch);
            DS2_info.DS2_peak_to_trough_amplitude(it_ch) = peak_to_trough_amplitude(it_ch);
            DS2_info.DS2_slope(it_ch) = peak_to_trough_slope(it_ch);
            if amplitude(it_ch) < -0.5 
                DS2_info.DS2_labels(it_ch) = "above";
            elseif amplitude(it_ch) >= -0.5 && amplitude(it_ch) <= 0.5
                DS2_info.DS2_labels(it_ch) = "in";
            elseif amplitude(it_ch) > 0.5
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
        load('single_shank_same_probe_map.mat','single_shank_same_probe_map')
        for it_ch = 1:32
            DS2_info.channel_depth(it_ch) = single_shank_same_probe_map(it_ch,3);
            %make amplitude and label features from selected spike 
            DS2_info.DS2_max_amplitude(it_ch) = amplitude(it_ch);
            DS2_info.DS2_mean_amplitude(it_ch) = mean_amplitude(it_ch);
            DS2_info.DS2_peak_to_trough_amplitude(it_ch) = peak_to_trough_amplitude(it_ch);
            DS2_info.DS2_slope(it_ch) = peak_to_trough_slope(it_ch);
            if amplitude(it_ch) < -0.5 
                DS2_info.DS2_labels(it_ch) = "above";
            elseif amplitude(it_ch) >= -0.5 && amplitude(it_ch) <= 0.5
                DS2_info.DS2_labels(it_ch) = "in";
            elseif amplitude(it_ch) > 0.5
                DS2_info.DS2_labels(it_ch) = "bellow";
            end
        end
    end

    
    save ([writeDir '\' char(dataset) '_DS2_info.mat'], 'DS2_info')
    

end 