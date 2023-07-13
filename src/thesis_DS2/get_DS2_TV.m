function get_DS2_TV( rat_ID, sleep_data,spike_idx, probe_type) 
% get_DS2 takes eegs from all chanels and filters them to detect DS2 spikes
% it then compares DS2 spike amplitude across channels to detect an
% inversion point (0 amp) and labels all channels at/around the point and
% above GCL and all below HL in a table 
% TO DO:
    % 1) Needs to make a matrix with DS2 labels per channel for any sleep 
    % trial thats loaded - columns : channel numbers that match spatData 
    % tetrodes (add maxampch identity to spatData in make wf features 
    % code), channel depth, DS2 voltages per channel (amplitude), DS2 label
    % per channel ("above", "in", "below"), distance from DS2 inversion in
    % absolute depth (only possible on shanks with inversion) 
    %2) find the sleep trial using info in spatData - add trial labels to the
    %table 
    %3) for now I have to be on the folder with the correct sleep file
    %that should change in the future. 

   
    % get sleep trial eegs
    sleeptrialname = sleep_data; 
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
    voltage_trace = eegs_S.*(1.5 ./ 128);
    voltage_trace = voltage_trace'; 
    
    %run extract spikes TW but put in Senzai and Buzsaki parameters for DS2
    %spikes 
    sample_rate = 250;
    spikeThreshold = 0.8;%1.14 based on Senzai and Buzsaki Threshold - this can be adjusted 
    

    % extract DS2 from eegs 
    [spike_mat,~]=extract_DS2_TV(voltage_trace,sample_rate,spikeThreshold);
    %[~, max_spread_idx] = max(spike_mat(1, 2, :),[],'all', 'linear');
    [~, spike_idx_sorted_by_spread] = sort(spike_mat(1, 2, :), 'descend');
    spike_idx = spike_idx_sorted_by_spread(spike_idx); %slect a putative ds2 spike to plot
    if probe_type == 1
        probe_data = load('single_shank_probe_map.mat');
    elseif probe_type == 2
        probe_data = load('multi_shank_probe_map.mat');
    elseif probe_type == 3
        probe_data = load('single_shank_same_probe_map.mat');
    end
    [probe_idx, probe_depth, unique_depths] = process_probe_data(probe_data);
    if probe_type == 1 || probe_type == 3 % single shank probe
        [amplitude_single, fig1, fig2] = plot_single_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat);
    elseif probe_type == 2 % multi shank probe
        [fig1]=plot_multi_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat);        
        [amplitude, mean_amplitude, fig2] = plot_ds2_inversion_line_plot(probe_depth, spike_idx, probe_idx, spike_mat);
    end
    
    %produce DS2 labels (in a table) per channel per rat from depth and voltage info
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
    
    sleep_trial = extractBefore(sleeptrialname, '.eeg');
    save (['C:\Users\Isabella\Documents\MATLAB\MatLabData\' rat_ID '_' sleep_trial '_DS2_info.mat'], 'DS2_info','fig1','fig2')

end

function [probe_idx, probe_depth, unique_depths] = process_probe_data(probe_data)
    % This function will preprocess probe data for the single shank case
    probe_data = struct2array(probe_data);
    probe_idx = probe_data(:, 1) + 1;
    probe_depth = probe_data(:, 3);
    unique_depths = length(unique(probe_depth));
end

function [amplitude_single, fig1, fig2]= plot_single_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat)
    % Plotting code for single shank.
    for k = 1: unique_depths
        cmap = colormap(bone(unique_depths));
        probe_color = probe_depth;
    end
    fig1 = figure();
    hold all
    cmap = colormap(bone(32));
    for i = 1:length(probe_idx)
       p_idx = probe_idx(i);
       a = plot(spike_mat(p_idx,3:end,spike_idx));
       a.Color = cmap(i, :);
    end

    % line plot 
    selected_spike = nan(2, 32);
    for i = 1:length(probe_idx)
        p_idx = probe_idx(i);
        selected_spike(1, i) = spike_mat(p_idx,19,spike_idx);
        selected_spike(2, i) = probe_depth(i);
    end

    fig2 = figure();
    plot(selected_spike(1,:),selected_spike(2,:));

    amplitude_single = selected_spike(1,:);
end

function [fig1] = plot_multi_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat)
    % Function to plot the multi shank. 4 figures, each 
    % one with a unique colour per depth (and a key). The x axis shows
    % voltage and the y axis shows time.
    for k = 1: unique_depths
        cmap = colormap(bone(unique_depths));
        probe_color = probe_depth;
    end
    fig1 = figure();
    for j = 1:4
        % This gives us the indices of each shank.
        subplot(2, 2, j)
        hold all
        for i = 1:8
           p_idx = probe_idx((j-1)*8+i);
           a = plot(spike_mat(p_idx,3:end,spike_idx));
           a.Color = cmap(i, :);
        end
        depth_str = cell(1, unique_depths);
        for i = 1:unique_depths
            depth_str{i} = num2str(probe_depth(i), '%.1f');
        end
        title("Plot of shank " + string(j))
        ylim([-1, 1])
        legend(depth_str);
    end
end

function [amplitude, mean_amplitude, fig2]= plot_ds2_inversion_line_plot(probe_depth, spike_idx, probe_idx, spike_mat)
    % Voltages and depths plotted per depth. Depth is on
    % the y-axis and voltage is on the x-axis. For the multi-shank case we
    % plot 4 subplots - one per shank.
    % line plot 
    amplitude =[]; %stitching the voltages back together into a 1x32 array for DS2_info table 
    mean_amplitude = []; %same for mean voltages 
    fig2 = figure();
    for j = 1:4
        % loop over shanks
        subplot(2, 2, j)
        hold all
        % loop over channels
        start_index = (j-1)*8 + 1;
        end_index = start_index + 7;
        p_indices = probe_idx(start_index:end_index);
        depths = probe_depth(1:8);
        voltages = spike_mat(p_indices,19,spike_idx);
        peak_amplitude_mean = mean(spike_mat(p_indices,17:21,spike_idx),2);
        plot(voltages, depths);
        title("Plot of voltage against depth for shank " + string(j))
        xlim([-1.5, 1.5])
        xlabel("Voltage")
        ylabel("Depth")
        amplitude = [amplitude;voltages];
        mean_amplitude = [mean_amplitude; peak_amplitude_mean];
    end
end
