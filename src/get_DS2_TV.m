function get_DS2_TV(spatData, spike_idx, probe_type)
% get_DS2 takes eegs from all chanels and filters them to detect DS2 spikes
% it then compares DS2 spike amplitude across channels to detect an
% inversion point (0 amp) and labels all channels at/around the point and
% above GCL and all below HL 

    %load all the eegs and put into a nsamp x 32 matrix called eegs - hc to 32
   
    
    %sleep trial eegs 
    
    %find the sleep trial using info in spatData - add trial labels to the
    %table
%     trialName = spatData.trialName;   
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
    spikeThreshold = 1.14;%based on Senzai and Buzsaki Threshold - this can be adjusted 
    

    % extract DS2 from eegs 
    [spike_mat,spike_count]=extract_DS2_TV(voltages,sample_rate,spikeThreshold);
    %[~, max_spread_idx] = max(spike_mat(1, 2, :),[],'all', 'linear');
    [~, spike_idx_sorted_by_spread] = sort(spike_mat(1, 2, :), 'descend');
    spike_idx = spike_idx_sorted_by_spread(spike_idx); %slect a putative ds2 spike to plot
    if probe_type == 1
        probe_data = load('single_shank_probe_map.mat');
    elseif probe_type == 2
        probe_data = load('multi_shank_probe_map.mat');
    end
    [probe_idx, probe_depth, unique_depths] = process_probe_data(probe_data);
    if probe_type == 1 % single shank probe
        plot_single_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat)

    elseif probe_type == 2 % multi shank probe
        plot_multi_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat)
        plot_ds2_inversion_line_plot(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat)
    end
    % we want to get the max, per spike window and its t_max
    % we want to get the value for t_max+1
end

function [probe_idx, probe_depth, unique_depths] = process_probe_data(probe_data)
    % This function will preprocess probe data for the single shank case
    probe_data = struct2array(probe_data);
    probe_idx = probe_data(:, 1) + 1;
    probe_depth = probe_data(:, 3);
    unique_depths = length(unique(probe_depth));
end

function plot_single_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat)
    % Plotting code for single shank.
    for k = 1: unique_depths
        cmap = colormap(bone(unique_depths));
        probe_color = probe_depth;
    end
    figure()
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

    figure()
    plot(selected_spike(1,:),selected_spike(2,:));
end

function plot_multi_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat)
    % Function to plot the multi shank. We want to create 4 figures, each 
    % one with a unique colour per depth (and a key). The x axis shows
    % voltage and the y axis shows time.
    for k = 1: unique_depths
        cmap = colormap(bone(unique_depths));
        probe_color = probe_depth;
    end
    figure;
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
        legend(depth_str);
    end
end

function plot_ds2_inversion_line_plot(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat)
    % We want to see the voltages and depths plotted per depth. Depth is on
    % the y-axis and voltage is on the x-axis. For the multi-shank case we
    % plot 4 subplots - one per shank.
    % line plot 
    figure()
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
        plot(voltages, depths);
        title("Plot of voltage against depth for shank " + string(j))
        xlim([-1.5, 1.5])
        xlabel("Voltage")
        ylabel("Depth")
    end
end

