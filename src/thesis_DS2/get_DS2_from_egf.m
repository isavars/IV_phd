function [max_amplitude, mean_amplitude, fig1, fig2, ds2_rate, peak_to_trough_amplitude, peak_to_trough_slope, spike_mat, new_eeg_chans]= get_DS2_from_egf(egf_data, probe_type, sws_table) 
% get_DS2 takes eegs from all chanels and filters them to detect DS2 spikes
% it also generates ds2 amplitude, ds2 rate, figures etc to analyze more
% features of ds2

%TO DO - 1) there's some hc in the amplitude making process - change this
%3) make dynamic artifact rejection - changes per age - troubleshoot by
%looking at all DS2s from r1311 2) make an if statement so that if
%sws_table is missing this can still run on an eeg file - you need to be
%able to make DS2's for days where you don't want to make a full cut 

%INPUTS: 1) egf_data - path to first egf in sleep trial,  2) probe_type - 2 
% is multishank, 1 and 3 are single shanks plugged in opposite ways and 4 
% is for tetrodes 3) sws_table - SWS from sleep
%filtering done on all the relevant trials. 
% EXAMPLE: get_DS2_from_egf('S:\DBIO_TKFC_SKGTIVA\thesis_data\PS_exposure_1\r1311\230213_1\230213h_sleepHP.egf', 2, 'C:\Users\Isabella\Documents\MATLAB\MatLabData\PS_1_exposure_sleepData.mat');

%OUTPUTS: useful features for DS2 analysis. 


    %parameters for voltage conversion and spike extraction 

    sample_rate = 4800; %in Hz
    spikeThreshold = 0.5;%1.14 Senzai and Buzsaki Threshold - this can be adjusted
    spikes_to_plot = 1; %this is just the best spike to make nice plots
%     best_DS2s = 50; %number of DS2s in spike_mat you want to include in DS2 info summary table (can change to all by becoming size(spike_mat,3) and being placed after spike_mat is made
    
    %load SWS data from table 
    load (sws_table, 'sleepData'); 

    %get trial name 
    lastBackslashIdx = find(egf_data == '\', 1, 'last');
    sleep_eeg = egf_data(lastBackslashIdx+1:end);
    sleep_trial = extractBefore(sleep_eeg, '.egf');
    %get current trial 
    curr_trial = find(sleepData.trialname == sleep_trial);
    %get dataset name 
    dataset = sleepData.dataset(curr_trial);
    %get rat age from dataset name 
    rat_age = str2double(extractAfter(dataset, 'P'));
    rat_id = str2double(extractBetween(dataset, 'r', '_'));
    

    %find number of files containing the current sleep trial name - this
    %makes it adaptabel for variable amounts of egf files 
    filePath = extractBefore(egf_data,sleep_trial);
    fileList = dir(fullfile(filePath, '*.egf*'));
    eeg_count = 0;
    eeg_ext = [];
    % Iterate through each file in the folder
    for i = 1:numel(fileList)
        % Extract the file name
        [~, fileName, fileExt] = fileparts(fileList(i).name);        
        % Check if the file name contains the search string
        if contains(string([fileName '.egf']), sleep_eeg)
            eeg_count = eeg_count + 1;
            file_ext = str2double(extractAfter( fileExt, '.egf'));
            if isnan(file_ext)
                eeg_ext = [eeg_ext; 1];
            else 
                eeg_ext = [eeg_ext; file_ext];
            end
        end
    end
   eeg_ext = sort(eeg_ext); %need to get in ascending order because the files are read in differently in fileList. 
    % get sleep trial eegs 
    sleeptrialname = egf_data;
    eegs_S = [];
    for it_eegs = 1:eeg_count
        if it_eegs == 1
            eeg_struct  = load_eeg(sleeptrialname);
        else
            eeg_struct = load_eeg(strcat(sleeptrialname,num2str(eeg_ext(it_eegs))));
        end
        eegs_S = [eeg_struct.eeg;eegs_S];
    end        
    eegs_S = reshape(eegs_S,[],eeg_count); 
    
    %r1242 has a noisy channel thats messing up the spike detection -
    %'ground' here. 
    if rat_id == 1242 
        eegs_S(:,23) = zeros(length(eegs_S(:,1)),1);
    else 
    end
    %convert to mV 
    voltage_trace = eegs_S.*(1.5 ./ 2^15); % divided by bit reolution x voltage on scope 

    %cut down sleep trials to sleep only  
    %find the sws epoch inds from the correct sleep file in the table 
    %sleepData.
    SWS_epochs = sleepData.SWS_epochs{curr_trial,1};
    SWS_epochs = SWS_epochs*sample_rate; %its possible that im adding samples in this conversion. 
    new_voltage_trace = [];
    for itEps = 1: length(SWS_epochs)
        start_index = int64(SWS_epochs(1,itEps));
        end_index = int64(SWS_epochs(2,itEps));
        new_voltage_trace = [new_voltage_trace; voltage_trace(start_index:end_index,:)];
    end 
    voltage_trace = new_voltage_trace; 

    %filter the data - two good filtering options made by TW 
    filtOrder     = 1000;
    filtObj       = designfilt( 'highpassfir','FilterOrder',filtOrder, 'CutoffFrequency', 10, 'SampleRate', 4800);
%         filtObj       = designfilt( 'highpassfir','FilterOrder',filtOrder, 'PassbandFrequency', 20, 'StopbandFrequency', 10, 'SampleRate', 4800);       
    voltage_trace = fftfilt( filtObj, voltage_trace );  % FFTFILT filters down *columns* of a 2D array
    voltage_trace = voltage_trace( (filtOrder/2):end, : );   % Need to correct for filter lag - remove order/2 samples from beginning of signal
    voltage_trace = [voltage_trace; nan(filtOrder/2, size(voltage_trace,2))];  % Pad with NaN to keep consistent size
    voltage_trace = voltage_trace';   % Now flip dims so matching expectations of rest of function.
   
    % extract DS2 from eegs 
    [spike_mat,~]=extract_DS2_from_egf(voltage_trace,sample_rate,spikeThreshold, rat_age); %changed to egf 
      
    %introduce a loop here that goes over best DS2s and makes means of
    %values to output for the bestDSs - it needs to break if spike_mat is
    %empty and it needs to make means of length spikemat 3rd dimenstion
    %instead of bestDS2s if its shorter than best DS2s 
    if ~isempty(spike_mat) %if extract_DS2 doesnt find any spikes that trial is skipped 

        %create spread sorted indexes for spike_mat 
        [~, spike_idx_sorted_by_spread] = sort(spike_mat(1, 2, :), 'descend'); 
        %make best_DS2s examples to be averaged in amplitude measures 
        num_spikes = size(spike_mat,3); %3rd dimension of the spikematrix is the number of spikes
        if num_spikes >= 100
            best_DS2s = round(num_spikes*0.2);% top 10 percentile of spikes will be the best quality - proportionate to amount detected 
        elseif num_spikes < 100 && num_spikes >= 20
            best_DS2s = 20;
        elseif num_spikes < 20
            best_DS2s = num_spikes;
        end 
        %make ds2_rate for analysis - numspikes/trial duration in seconds 
        trial_duration = sleepData.duration(curr_trial);
        ds2_rate = num_spikes/trial_duration; 
        %load probe maps and probe data to be used in figure making functions 
        if probe_type == 1
            probe_data = load('single_shank_probe_map.mat');
            new_eeg_chans = [];
        elseif probe_type == 2
            probe_data = load('multi_shank_probe_map.mat');
            new_eeg_chans = [];
        elseif probe_type == 3
            probe_data = load('single_shank_same_probe_map_final.mat'); 
            new_eeg_chans = [];
        elseif probe_type == 4 
            %load eeg_chans and chop down spike_mat if the data is not from probes 
            eeg_chans = sleepData.eeg_channels(curr_trial); %maybe can use these to color the plots     
            [new_eeg_chans, spike_mat] = remove_channel_repeats(eeg_chans, spike_mat);
            probe_data = [];
        end
        if ~isempty(probe_data)
            [probe_idx, probe_depth, unique_depths] = process_probe_data(probe_data);
        end
        % make max_amplitude, mean_amplitude and figures
        if size(spike_mat,3) >= spikes_to_plot       
            spike_idx_for_best_spikes = [];
            for ii = 1:best_DS2s
                spike_idx_for_best_spikes = [spike_idx_for_best_spikes, spike_idx_sorted_by_spread(ii)];
            end 
             %select best ds2 spikes to plot             
            for jj = 1:spikes_to_plot
                spike_idx = spike_idx_sorted_by_spread(jj);                
                %call figure making functions and feature making functions       
                %figure making functions: if you want to plot more than one this
                %should be in a loop length spikes_to_plot 
                if probe_type == 1 || probe_type == 3 % single shank probe
                    [fig1, fig2] = plot_single_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat, dataset);
                elseif probe_type == 2 % multi shank probe
                    [fig1] = plot_multi_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat, dataset);        
                    [fig2] = plot_ds2_inversion_line_plot(probe_depth, spike_idx, probe_idx, spike_mat);
                elseif probe_type == 4 %tetrodes (repeat below)
                    %create a plot of the DS2s kind of like the single shank
                    %plot but with no mapping
                    [fig1] = plot_tetrodes(spike_mat, spike_idx, dataset, new_eeg_chans);
                    fig2= nan;
                 end
             end 

            %call feature funciton 
            [max_amplitude, mean_amplitude, peak_to_trough_amplitude, peak_to_trough_slope]= make_DS2_features(spike_mat, spike_idx_for_best_spikes);

        else %make features and figures from average of availabe spikes if there are very few spikes
            spike_idx_for_best_spikes = [];
            for ii = 1:num_spikes
                spike_idx_for_best_spikes = [spike_idx_for_best_spikes, spike_idx_sorted_by_spread(ii)];
            end 
            
            for jj = 1:num_spikes %should be capped but its ok since its under 10
                spike_idx = spike_idx_sorted_by_spread(jj);
                if probe_type == 1 || probe_type == 3 % single shank probe
                    [fig1, fig2] = plot_single_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat, dataset);
                elseif probe_type == 2 % multi shank probe
                    [fig1]=plot_multi_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat, dataset);        
                    [fig2] = plot_ds2_inversion_line_plot(probe_depth, spike_idx, probe_idx, spike_mat);
                elseif probe_type == 4 %tetrodes (repeat below)
                    %create a plot of the DS2s kind of like the single shank
                    %plot but with no mapping  
                    [fig1] = plot_tetrodes(spike_mat, spike_idx, dataset, new_eeg_chans);
                    fig2= nan;
                end  
            end
            %call feature funciton 
            [max_amplitude, mean_amplitude, peak_to_trough_amplitude, peak_to_trough_slope]= make_DS2_features(spike_mat, spike_idx_for_best_spikes);
        end 

   else  %produce empty values to fill that DS2_info 
        max_amplitude = nan(32,1);
        mean_amplitude = nan(32,1);
        ds2_rate = nan;
        peak_to_trough_amplitude = nan(32,1);
        peak_to_trough_slope = nan(32,1);
        fig1 = nan;fig2= nan; %this seems like it should be done some other way
        eeg_chans = sleepData.eeg_channels(curr_trial); %maybe can use these to color the plots     
        [new_eeg_chans] = remove_channel_repeats(eeg_chans);
    end

end

function [probe_idx, probe_depth, unique_depths] = process_probe_data(probe_data)
    % This function will preprocess probe data for the single shank case
    probe_data = struct2array(probe_data);
    probe_idx = probe_data(:, 1) + 1;
    probe_depth = probe_data(:, 3);
    unique_depths = length(unique(probe_depth));
end

function [fig1, fig2]= plot_single_shank(probe_depth, unique_depths, spike_idx, probe_idx, spike_mat,dataset)
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
    title(dataset)

    % line plot 
    selected_spike = nan(3, 32);
    for i = 1:length(probe_idx)
        p_idx = probe_idx(i);
        selected_spike(1, i) = spike_mat(p_idx,339,spike_idx);
        selected_spike(2, i) = probe_depth(i);
        selected_spike(3, i) = mean(spike_mat(p_idx,321:353,spike_idx));
    end

    fig2 = figure();
    plot(selected_spike(1,:),selected_spike(2,:));
    title(dataset)

end

function [fig1] = plot_tetrodes(spike_mat, spike_idx, dataset, new_eeg_chans)
    
    %make colors in the theme 
    cmap = colormap(parula(length(new_eeg_chans)));
    %make legend values match the tetrode instead of the channel
    chan_tets = zeros(64, 2);
    chan_tets(:, 1) = (1:64)';
    chan_tets(:, 2) = repelem(1:16, 4)';
    % Find the matching indices in the first column of chan_tets
    [~, indices] = ismember(new_eeg_chans, chan_tets(:, 1));
    % Replace the values in new_eeg_chans with the corresponding values from the second column of chan_tets
    new_eeg_chans_tet = chan_tets(indices, 2);

    %make plot 
    fig1 = figure();
    hold all
    for i = 1:length(new_eeg_chans)
        a = plot(spike_mat(i,3:end,spike_idx));
        a.Color = cmap(i, :);
    end
    title(dataset)
    legend(num2str(new_eeg_chans_tet))
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
        ylim([-1.5, 1.5])
        xlim([0 700])
        legend(depth_str);
    end
end

function [fig2]= plot_ds2_inversion_line_plot(probe_depth, spike_idx, probe_idx, spike_mat)
    % Voltages and depths plotted per depth. Depth is on
    % the y-axis and voltage is on the x-axis. For the multi-shank case we
    % plot 4 subplots - one per shank.
    %it also makes the amplitude measures to be used in analysis further
    %along

    fig2 = figure();
    for j = 1:4
        % loop over shanks
        subplot(2, 2, j);
        hold all
        % loop over channels
        start_index = (j-1)*8 + 1;
        end_index = start_index + 7;
        p_indices = probe_idx(start_index:end_index);
        depths = probe_depth(1:8);
        voltages = spike_mat(p_indices,339,spike_idx); %changes to hc values for egf sampling rate (this is the middle sample)
        plot(voltages, depths);
        title("Plot of voltage against depth for shank " + string(j))
        xlim([-1.5, 1.5])
        xlabel("Voltage")
        ylabel("Depth")

    end
end

function [max_amplitude, mean_amplitude, peak_to_trough_amplitude, peak_to_trough_slope]= make_DS2_features(spike_mat, spike_idx_for_best_spikes)
        %make max_ amplitude, mean_amplitude, peak_to_trough amplitude and
        %slope means for x many DS2s found - need to be 1x32 for DS2_info

        %make num_eegs so its compatible with different amounts of eegs 
        num_eegs = size(spike_mat,1); %check this is good 
        
        %loop over best DS2s and produce mean values of all the features
        max_amplitudes = [];
        mean_amplitudes = [];
        peak_to_trough_amplitudes = [];
        peak_to_trough_slopes = [];

        for spike_num = spike_idx_for_best_spikes
            % get the max point in the voltage trace for all the channels then
            % search the next 10ms for a trough
            peaks = spike_mat(:,339, spike_num); % matrix of voltage values at spike peak (centre of the window)
            second_half_of_spike = spike_mat(:,339:339+48, spike_num);
            [~,max_amp_chan] = max(abs(peaks));
            [~, trough_idx] =  min(abs(second_half_of_spike(max_amp_chan,:)),[],2);
            troughs = second_half_of_spike(:, trough_idx);
            peak_to_trough_amplitude = peaks - troughs; %incase the spike is riding on another spike
            peak_to_trough_time = trough_idx ./ 4800; %convert to seconds 
            peak_to_trough_slope = (troughs - peaks) ./ peak_to_trough_time; %calculate slope from peak to trough -best for deteriming spike orientation
            peak_amplitude_mean = mean(spike_mat(:,332:342,spike_num),2); %this corresponds to the middle 10% of the spike - trying to smooth the peak here
      
            %accumulate spike info into a matrix 
            max_amplitudes = [max_amplitudes, peaks]; %peaks is amplitude at max value
            mean_amplitudes = [mean_amplitudes, peak_amplitude_mean];
            peak_to_trough_amplitudes = [peak_to_trough_amplitudes, peak_to_trough_amplitude];
            peak_to_trough_slopes = [peak_to_trough_slopes, peak_to_trough_slope]; 
        end
        max_amplitude = mean(max_amplitudes,2);
        mean_amplitude = mean(mean_amplitudes,2);
        peak_to_trough_amplitude = mean(peak_to_trough_amplitudes,2);
        peak_to_trough_slope = mean(peak_to_trough_slopes,2);

        if num_eegs < 32 
            empty_chans = nan((32 - num_eegs),1);
            max_amplitude = [max_amplitude; empty_chans];
            mean_amplitude = [mean_amplitude; empty_chans];
            peak_to_trough_amplitude = [peak_to_trough_amplitude; empty_chans];
            peak_to_trough_slope = [peak_to_trough_slope; empty_chans];
        else 
        end 
end 
function [new_eeg_chans, new_spike_mat] = remove_channel_repeats(eeg_chans, spike_mat)   
    % remove repearted channels - if the same eeg is found twice plot the 
    % last one because the ones made on direct should be in the first 
    % positions and when eeg1 is one of the channel it will grab the next 
    % one but if tis the last one its sorted 
    
    %condition if spike_mat is not there
    if nargin < 2
        spike_mat = [];
    end
    %convert eeg_chans into a numerical array 
    eeg_chans = cell2mat(eeg_chans{1});
    % Find indices of repeated values
    [~, ~, idx] = unique(eeg_chans, 'stable');    
    % Find indices of the last occurrence of each unique value
    last_indices = find(histcounts(idx) > 1);    
    % Create a logical array to mark only the last occurrence of each value
    is_last_occurrence = true(size(eeg_chans));
    is_last_occurrence(last_indices) = false;
    new_eeg_chans = eeg_chans(is_last_occurrence);    
    if ~isempty(spike_mat)
        new_spike_mat = spike_mat(is_last_occurrence,:,:);
    else
        new_spike_mat = [];
    end

end