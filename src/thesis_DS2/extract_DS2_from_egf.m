function [spike_mat,spike_count] = extract_DS2_from_egf(eeg_data,samp_rate,spikeThreshold, rat_age)
% Extract spike data (150ms window) from all channels and extract corresponding
% timestamp. Spikes are considered voltages which cross the DACQ threshold
% sample rate for egf is 4800Hz so 140ms is 672 samples - the spikes (~20ms)
% should occur over 96 samples. 

% INPUT:
% egf_data: number_of_egfs*samples matrix of voltage data from all channels in sleep
% trial
% samp_rate: sampling rate 4800 Hz
% spikeThreshold: going with 0.8 for now. 
% OUTPUT:
% spike_count: total number of dentate spikes recorded
% spike_mat: matrix containing spike time (row 1) and corresponding
%   voltages across all channels.
% TO DO: 1) make adaptable artifact rejection based on age of the rat -
% will need to play around with this a bit. 2) there needs to be another
% input to this function - rat age - can get it from sleepData.dataset 

eeg_data_abs= abs(eeg_data);%get the absolute value so positive and negative threshold values are considered. 

thr_cr         = max(eeg_data_abs,[],1) > spikeThreshold;  % 'thr_cr' is a logical index of where threshold is crossed
thr_cr_numind  = find( thr_cr );                  % 'thr_cr_numind' is a numerical index for where threshold crossed

possible_count    = numel(thr_cr_numind);              % Maximum possible no. spikes. 
nSamp             = size( eeg_data, 2 );

%set hc values for spike regions of interest
trace_window = [-336 336]; %set trace window of ~140ms (673 samples so its an odd number for centering)
spike_window = 289:384;  % Specify the column range of interest - 20ms in the center of the window 

spike_count = 0;
thr_cr_idx  = find( thr_cr_numind>trace_window(2), 1, 'first' );    % 'thr_cr_idx' is the iterator over 'thr_cr_numind' for the loop below. Get the earliest >sample 336, so as to allow pre-crossing samples.

% Loop through crossings and extract spikes

spike_mat = cell(1, 3); %cell array that grows helps to deal with the larger file size of the egf. 
    
while thr_cr_idx <= length(thr_cr_numind) 
 
    thr_cr_samp = thr_cr_numind(  thr_cr_idx  );  % 'thr_cr_samp' is the current number of the sample where threshold is crossed.
    
    if thr_cr_samp > nSamp-trace_window(2) || thr_cr_samp <= trace_window(2)  % If we have reached the end of the data or we are too close to the beginning
        break
    else 

        put_spk_data = eeg_data_abs( :, thr_cr_samp + (trace_window(1):trace_window(2)) );  % 'put_spk_data' = putative spike data, putative as in hasn't undergone artifact filtering yet.
        raw_put_spk_data = eeg_data( :, thr_cr_samp + (trace_window(1):trace_window(2)) );
        % Get the necessary info for testing artifact rejection
        
        [max_value, max_ind_in_columns] = max(put_spk_data(:, spike_window), [], 'linear');
        max_val= max(max_value);
        max_val_column_idx = find (max_value == max_val,1, 'first');
        max_ch = max_ind_in_columns(1);
        max_timestamp_idx = spike_window(max_val_column_idx);

        % make variables for artifact filtering - not sure if these are
        % centered around the peak of the spike...
        first_samps        = put_spk_data( max_ch, spike_window(1:length(spike_window)/3)); %this is the first samples in the spike (1st 3rd of the spike)
        last_samps       = put_spk_data( max_ch, spike_window(end-length(spike_window)/3:end));% last samples in the spike
        
        % try to put something in that finds the dip after the spike? 
        
        %make spread array 
        spread_array = raw_put_spk_data(:, max_timestamp_idx);
        spread = max(spread_array) - min(spread_array);

        % we believe that this array represents a single trace for a
        % putative spike.  %check the spike has the right shape
        putative_spike_window = put_spk_data(max_ch, :);
        [passes_waveform_shape_filter] = waveform_shape_filter(putative_spike_window, max_timestamp_idx, max_val);
                              
        % needs to remove things that are not spike shaped from the section that should contain the spike 
        % Check if all the conditions are met

        %add if statement for rat age to create a couple different filters
        %here:
        if rat_age >= 22
            age_based_spread_cutoff = 0.7;
        elseif rat_age < 22 && rat_age >= 20
            age_based_spread_cutoff = 0.5;
        elseif rat_age < 20
            age_based_spread_cutoff = 0.2;
        end 
        age_based_spread_cutoff = 0.2; %to let as many in to not add an artificial developmental effect (for amplitude ds2 plot)

        if  passes_waveform_shape_filter && mean(first_samps) < max_val*0.30  &&  mean(last_samps) < max_val*0.30 && spread > age_based_spread_cutoff % This is artefact rejection test - can modify these numbers            
            % Assume it is a real spike and store
            spike_count = spike_count + 1;           
            max_voltage_idx = thr_cr_samp + max_timestamp_idx - trace_window(2);
            max_put_spk_data = eeg_data(:, max_voltage_idx + (trace_window(1):trace_window(2))); % now that we have the max index between timesteps of the threshold we redefine put_spk_data to be centered on the max
            spike_mat{:,1,spike_count}     = max_voltage_idx/samp_rate*10^6; % Save spike timestamp for all channels - in microseconds - added by IV
            spike_mat{:,2,spike_count}     = spread; % Save spread of the spike
            spike_mat{:,3:end,spike_count} = max_put_spk_data; %trying something - curly brakets on all 
        end

        thr_cr_idx  = find( thr_cr_numind>(thr_cr_samp+673), 1, 'first' );

        if isempty( thr_cr_idx )
            break
        end

    end

end

%sow this back into a matrix so it works in the code that follows 
num_eegs = size(eeg_data,1);% get the number of eegs recorded for the trial
matrix_spike_mat = zeros(num_eegs, 675, spike_count); 

for i = 1:spike_count
    timestamp = spike_mat{:,1,i};
    spread = spike_mat{:,2,i};
    eeg_sample = spike_mat{:,3,i};
    
    matrix_spike_mat(:, 1, i) = repmat(timestamp, num_eegs, 1);
    matrix_spike_mat(:, 2, i) = repmat(spread, num_eegs, 1);
    matrix_spike_mat(:, 3:end, i) = eeg_sample;
end

spike_mat = matrix_spike_mat;
whos
end

function [passes_waveform_shape_filter] = waveform_shape_filter(trace, max_timestamp_idx,max_val)
    % This function takes in an array and returns a logical
    spike_window = trace(max_timestamp_idx - 48: max_timestamp_idx +48); % 24 samples at either side of a the spike sample window which gives a tolerance +/- 5ms to the width to determine if its too wide (more than 30ms)
%     baseline_range = mean(put_spk_data([1:226 408:end])); % mean of the voltage trace at either side of the max spike window 
%     min_spike_width = 24; % 5ms in samples 
%     start_index = ceil((length(max_spike_window) - min_spike_width) / 2) +1 ;
%     end_index = start_index + min_spike_width -1;
%     middle_indices_max_spike_window = start_index : end_index;
%     min_spike_window = max_spike_window(middle_indices_max_spike_window); %less hc than the max - do this to the rest of the code 
%     middle_index_min_spike_window = min_spike_width/2 +0.5;
    values_over_half_height = spike_window > max_val*0.5;
    % Calculate the middle index of the array
    middleIndex = ceil(length(values_over_half_height) / 2);
    % Find the change points in the array
    changePoints = diff([false,values_over_half_height,false]);
    % Identify the start and end indices of each segment
    startIndices = find(changePoints==1);
    endIndices = find(changePoints==-1);
    % Initialize an empty array to hold the lengths of the true segments
    trueSegmentLengths = [];
    % Initialize an empty array to hold the start indices of the true segments
    trueSegmentStarts = [];
    % Loop through the segments
    for i = 1:length(startIndices)
        % If the segment is a true segment and it contains the middle index, 
        % save the start index and length of the segment
        if startIndices(i) <= middleIndex && endIndices(i) >= middleIndex
            trueSegmentStarts = [trueSegmentStarts; startIndices(i)];
            trueSegmentLengths = [trueSegmentLengths; endIndices(i) - startIndices(i)];
        end
    end
    if ~isempty(trueSegmentLengths)
        half_height_width = trueSegmentLengths; 
        passes_waveform_shape_filter =  half_height_width < 72  && half_height_width > 12; %half height width should be half the expected width at the base  
        % Check if the spike width is all above baseline at the narrowest option and not a max number wider at either side of the expected 20ms so its excluding any spikes wider than 36ms            
    else 
        passes_waveform_shape_filter = false; 
    end
end 

