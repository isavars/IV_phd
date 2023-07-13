function [amplitude, mean_amplitude, fig1, fig2, ds2_rate]= get_DS2_IV(eeg_data,spike_num, probe_type, sws_table,curr_trial) 
% get_DS2 takes eegs from all chanels and filters them to detect DS2 spikes
% it also generates ds2 amplitude, ds2 rate, figures etc to analyze more
% features of ds2

% TO DO: 
    %1) amplitude needs to be a trough to peak measure not the voltage off
    %the oscilloscope. 
    %2) add egf insread of eeg (make sure all the egf files are joined also
    %and in the right folders) 
    %3) if spike extraction fails the trial needs to be skipped 

    %parameters for voltage conversion and spiken extraction 

    sample_rate = 250;
    spikeThreshold = 0.8;%1.14 based on Senzai and Buzsaki Threshold - this can be adjusted 
    
    %load SWS data from table 
    load (sws_table, 'sleepData'); 

    %get trial name 
    lastBackslashIdx = find(eeg_data == '\', 1, 'last');
    sleep_eeg = eeg_data(lastBackslashIdx+1:end);
    sleep_trial = extractBefore(sleep_eeg, '.eeg');
    %get current trial 
    curr_trial = find(sleepData.trialname == sleep_trial);
    %get dataset name 
    dataset = sleepData.dataset(curr_trial);
   
    % get sleep trial eegs 
    sleeptrialname = eeg_data;
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

    %cut down sleep trials to sleep only  
    %find the sws epoch inds from the correct sleep file in the table 
    %sleepData.
    SWS_epochs = sleepData.SWS_epochs{curr_trial,1};
    %convert epoch inds into whatever the voltage trace is in or actually
    %maybe it makese sense to add it as a filter inside of extract_DS2 once
    %its in seconds? 
    SWS_epochs = SWS_epochs*sample_rate; %its possible that im adding samples in this conversion. 
    new_voltage_trace = [];
    for itEps = 1: length(SWS_epochs)
        start_index = int64(SWS_epochs(1,itEps));
        end_index = int64(SWS_epochs(2,itEps));
        new_voltage_trace = [new_voltage_trace, voltage_trace(:,start_index:end_index)];
    end 
    voltage_trace = new_voltage_trace; 
   
    % extract DS2 from eegs 
    [spike_mat,~]=extract_DS2_IV(voltage_trace,sample_rate,spikeThreshold);
    if ~isempty(spike_mat) %if extract_DS2 dosnt find any spikes that trial is skipped
        [~, spike_idx_sorted_by_spread] = sort(spike_mat(1, 2, :), 'descend');
        spike_idx = spike_idx_sorted_by_spread(spike_num); %select a putative ds2 spike to plot
        if probe_type == 1
            probe_data = load('single_shank_probe_map.mat');
        elseif probe_type == 2
            probe_data = load('multi_shank_probe_map.mat');
        elseif probe_type == 3
            probe_data = load('single_shank_same_probe_map.mat');
        end
        [probe_idx, probe_depth, unique_depths] = process_probe_data(probe_data);
        if probe_type == 1 || probe_type == 3 % single shank probe
            [amplitude,mean_amplitude, fig1, fig2] = plot_single_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat, dataset);
        elseif probe_type == 2 % multi shank probe
            [fig1]=plot_multi_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat, dataset);        
            [amplitude, mean_amplitude, fig2] = plot_ds2_inversion_line_plot(probe_depth, spike_idx, probe_idx, spike_mat);
        end
        
        %make ds2_rate for analysis - numspikes/trial duration in seconds 
        num_spikes = size(spike_mat,3); %3rd dimension of the spikematrix is the number of spikes 
        trial_duration = sleepData.duration(curr_trial);
        ds2_rate = num_spikes/trial_duration; 

    else %produce empty values to fill that DS2_info 
        amplitude = nan(32,1);
        mean_amplitude = nan(32,1);
        ds2_rate = nan;
        fig1 = nan;fig2= nan; %this seems like it should be done some other way
    end
end

function [probe_idx, probe_depth, unique_depths] = process_probe_data(probe_data)
    % This function will preprocess probe data for the single shank case
    probe_data = struct2array(probe_data);
    probe_idx = probe_data(:, 1) + 1;
    probe_depth = probe_data(:, 3);
    unique_depths = length(unique(probe_depth));
end

function [amplitude, mean_amplitude, fig1, fig2]= plot_single_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat,dataset)
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
    selected_spike = nan(3, 32);
    for i = 1:length(probe_idx)
        p_idx = probe_idx(i);
        selected_spike(1, i) = spike_mat(p_idx,19,spike_idx);
        selected_spike(2, i) = probe_depth(i);
        selected_spike(3, i) = mean(spike_mat(p_idx,17:21,spike_idx));
    end

    fig2 = figure();
    plot(selected_spike(1,:),selected_spike(2,:));
    title(dataset)

    amplitude = selected_spike(1,:)';
    mean_amplitude = selected_spike(3,:)';
end

function [fig1] = plot_multi_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat, dataset)
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
        title(dataset + " shank " + string(j))
        ylim([-1, 1])
        legend(depth_str);
    end
end

function [amplitude, mean_amplitude, fig2]= plot_ds2_inversion_line_plot(probe_depth, spike_idx, probe_idx, spike_mat)
    % Voltages and depths plotted per depth. Depth is on
    % the y-axis and voltage is on the x-axis. For the multi-shank case we
    % plot 4 subplots - one per shank.
    %it also makes the amplitude measures to be used in analysis further
    %along

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

      