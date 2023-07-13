function get_DS2_info_wrapper(num_spikes, sleep_data, probe_type, writeDir, sws_table)
%I want this function to loop through get_DS2_IV and collect the first n
%DS2s and give mean values to make DS2 info. Then makes and saves DS2 info for
%elePos to be used in classification and DS2 analysis in general. 
% Makes a matrix with DS2 labels per channel for any sleep 
    % trial thats loaded - columns : channel numbers that match spatData 
    % tetrodes, channel depth, DS2 voltages per channel (amplitude), DS2 label
    % per channel ("above", "in", "below"), distance from DS2 inversion in
    % absolute depth (only possible on shanks with inversion) 
%TO DO 1) add something that grabs all the files called DS2_info in a
%folder and makes electrode position tables from it so just transposes the
%data and adds it to a table but gets the rat ID and age per row from bits
%here. 2) reducde the number of input arguments - numspikes can become a
%parameter when you decide on a good one. the other things could be read
%off a spreadsheet. writter dir and sws_table should stay 3) add DS2 rate
%-can be picked up by get_DS2_IV

    %get trial name for saving files (probably there's a better way to do
    %this) 
    lastBackslashIdx = find(sleep_data == '\', 1, 'last');
    sleep_eeg = sleep_data(lastBackslashIdx+1:end);
    sleep_trial = extractBefore(sleep_eeg, '.eeg');
    %load SWS data from table to get parameters 
    load (sws_table, 'sleepData');  
    curr_trial = find(sleepData.trialname == sleep_trial);
    dataset = sleepData.dataset(curr_trial);

    %call get_DS2_IV to get all the necessary bits 
    amplitudes = [];
    mean_amplitudes = [];
    for spike_idx = 1:num_spikes
        [amplitude, mean_amplitude, fig1, fig2]= get_DS2_IV( sleep_data,spike_idx, probe_type,sws_table);
        amplitudes = [amplitudes, amplitude];
        mean_amplitudes = [mean_amplitudes, mean_amplitude];
        cur_spike = num2str(spike_idx);
        %savefig ([writeDir '\' rat_ID '_' sleep_trial '_DS2_spike_figure_' cur_spike '.mat'],'fig1','fig2')
    end 
    amplitude = mean(amplitudes,2);
    mean_amplitude = mean(mean_amplitudes,2);

    %produce DS2 labels (in a table) per channel per rat per day from depth and voltage info
    %used in the inversion line plots. makes 3 kinds of table depending on
    %probe type. 

    % Define the variable names and types
    varNames = {'channel', 'channel_depth', 'DS2_labels', 'DS2_amplitude', 'DS2_mean_amplitude'}; %'DS2_distance'};
    varTypes = {'double', 'double', 'string', 'double', 'double'};

    % Create an empty table with the specified shape
    DS2_info = table('Size', [32, 5], 'VariableNames', varNames, 'VariableTypes', varTypes);
    %add channel numbers to the table (same for all probes) 
    DS2_info.channel = (1:32)'; 
    
    % insert data into the table 
    if probe_type == 1 
        load('single_shank_probe_map.mat','single_shank_probe_map')
        DS2_info.channel_depth(it_ch) = single_shank_probe_map(it_ch,3);
        %make features from selected spike 
        DS2_info.DS2_amplitude(it_ch) = amplitude_single(it_ch); 
        if amplitude_single(it_ch) < -0.5 
            DS2_info.DS2_labels(it_ch) = "above";
        elseif amplitude_single(it_ch) >= -0.5 && amplitude_single(it_ch) <= 0.5
            DS2_info.DS2_labels(it_ch) = "in";
        elseif amplitude_single(it_ch) > 0.5
            DS2_info.DS2_labels(it_ch) = "bellow";
        end
        
    elseif probe_type == 2
        load('multi_shank_probe_map.mat','multi_shank_probe_map')
        for it_ch = 1:32              
            DS2_info.channel_depth(it_ch) = multi_shank_probe_map(it_ch,3); %check that mapping makes sense - does this need to be converted? 
            DS2_info.DS2_amplitude(it_ch) = amplitude(it_ch); % fill amplitudes per channel just by indexing the voltages 
            DS2_info.DS2_mean_amplitude(it_ch) = mean_amplitude(it_ch);
            % amplitude values chosen based on one good DS2 example - can
            % tweak in the future. 
            if amplitude(it_ch) < -0.5 
                DS2_info.DS2_labels(it_ch) = "above";
            elseif amplitude(it_ch) >= -0.5 && amplitude(it_ch) <= 0.5
                DS2_info.DS2_labels(it_ch) = "in";
            elseif amplitude(it_ch) > 0.5
                DS2_info.DS2_labels(it_ch) = "bellow";
            end
            % DS2_info.DS2_distance = leaving for another day - will
            % need to label the inflection point 0 um and then only label
            % chanels on shanks with inflection points with distance from
            % that 0 in um using 'channel_depth'
        end
    elseif probe_type ==3 
        load('single_shank_same_probe_map.mat','single_shank_same_probe_map')
        DS2_info.channel_depth(it_ch) = single_shank_same_probe_map(it_ch,3);
        %make amplitude and label features from selected spike 
        DS2_info.DS2_amplitude(it_ch) = amplitude(it_ch); 
        if amplitude(it_ch) < -0.5 
            DS2_info.DS2_labels(it_ch) = "above";
        elseif amplitude(it_ch) >= -0.5 && amplitude(it_ch) <= 0.5
            DS2_info.DS2_labels(it_ch) = "in";
        elseif amplitude(it_ch) > 0.5
            DS2_info.DS2_labels(it_ch) = "bellow";
        end
    end

        
    save ([writeDir '\' dataset '_DS2_info.mat'], 'DS2_info')


end 