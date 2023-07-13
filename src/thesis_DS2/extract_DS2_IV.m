function [spike_mat,spike_count] = extract_DS2_IV(eeg_data,samp_rate,spikeThreshold)
% Extract spike data (50ms window) from all channels and extract corresponding
% timestamp. Spikes are considered voltages which cross the DACQ threshold
% (1.14 V - is value from Senzai and Buzsaki - lowering for pups).
% 140ms - 35 samples
% INPUT:
% eeg_data: 32*samples matrix of voltage data from all channels in sleep
% trial
% OUTPUT:
% spike_count: total number of dentate spikes recorded
% spike_mat: matrix containing spike time (row 1) and corresponding
%   voltages across all channels.
% TO DO: 1) consider making this a subfunction of get_DS2 2) need to
% optimize filtering for spikes. 

eeg_data_abs= abs(eeg_data);%get the absolute value so positive and negative threshold values are considered. 

thr_cr         = max(eeg_data_abs,[],1) > spikeThreshold;  % 'thr_cr' is a logical index of where threshold is crossed
thr_cr_numind  = find( thr_cr );                  % 'thr_cr_numind' is a numerical index for where threshold crossed

possible_count    = sum(thr_cr);              % Maximum possible no. spikes. 
spike_mat         = nan(32,37,possible_count); % Matrix to store spike data
nSamp             = size( eeg_data, 2 );

% upperThreshold    = 250; % this is legacy - maybe there is no need for an upper threshold 

spike_count = 0;
thr_cr_idx  = find( thr_cr_numind>15, 1, 'first' );    % 'thr_cr_idx' is the iterator over ''thr_cr_numind' for the loop below. Get the earliest >sample 10, so as to allow pre-crossing samples.

% Loop through crossings and extract spikes
while thr_cr_idx <= length(thr_cr_numind)
    
    thr_cr_samp = thr_cr_numind(  thr_cr_idx  );  % 'thr_cr_samp' is the current number of the sample where threshold is crossed.
    
    if thr_cr_samp > nSamp-17 || thr_cr_samp <= 17  % If we have reached the end of the data or we are too close to the beginning
        break
    else 

        put_spk_data = eeg_data_abs( :, thr_cr_samp + (-17:17) );  % 'put_spk_data' = putative spike data, putative as in hasn't undergone artifact filtering yet.
        raw_put_spk_data = eeg_data( :, thr_cr_samp + (-17:17) );
        % Get the necessary info for testing artifact rejection
%         [max_val,max_ind] = max(put_spk_data,[],'all', 'linear');
%         [max_ch, max_timestamp_idx]= ind2sub( size(put_spk_data), max_ind);
        spike_index = 15:19;  % Specify the column range of interest
        [max_value, max_ind_in_columns] = max(put_spk_data(:, spike_index), [], 'linear');
        max_val= max(max_value);
        max_val_column_idx = find (max_value == max_val,1, 'first');
        max_ch = max_ind_in_columns(1);
        max_timestamp_idx = spike_index(max_val_column_idx);

        % make sections for artifact filtering 
        first_samps        = put_spk_data( max_ch, 14:15); %this is now the first samples in the spike not in the full 140ms window
        last_samps       = put_spk_data( max_ch, 19:20);% last samples in the spike
        % try to put something in that finds the dip after the spike? 
        %make spread array 
        spread_array = raw_put_spk_data(:, max_timestamp_idx);
        spread = max(spread_array) - min(spread_array);
 
       % Check spike waveform shape
        max_spike_window = put_spk_data(max_ch, 13:21); % 2 samples at either side of a 7 sample window which gives a tolerance +/- 8ms to the width to determine if its too wide
        baseline_range = mean(put_spk_data(max_ch, [1:13 21:end]));
        min_spike_width = 3;
        start_index = ceil((length(max_spike_window) - min_spike_width) / 2) +1 ;
        end_index = start_index + min_spike_width -1;
        middle_indices = start_index : end_index;
        min_spike_window = max_spike_window(middle_indices); %less hc than the max
        
        % Check if the spike width is a max of one sample wider on either side of the expected 20ms so its excluding any spikes wider than 36ms            
        is_within_width_range = all(min_spike_window > baseline_range); %&& ~any((max_spike_window(1) || max_spike_window(end)) > max_val*0.33) ;             

        % Check if the maximum amplitude (max_val) matches the desired amplitude within tolerance
        desired_amplitude = 1.5; % Update with the desired amplitude value
        amplitude_tolerance = 1; % Update with the amplitude tolerance value
        is_within_amplitude_tolerance = (max_val >= (desired_amplitude - amplitude_tolerance)) && (max_val <= (desired_amplitude + amplitude_tolerance)); % probably just has to be bigger than the min - is this the same as the threshold??
        % needs to remove things that are not spike shaped from the section that should contain the spike (two samples before and after the threshold crossing)
        if ~any(first_samps > max_val*0.33)  &&  ~any(last_samps > max_val*0.33) % This is artefact rejection test - can modify these numbers            
            % Check if all the conditions are met
%                 if  is_within_width_range && is_within_amplitude_tolerance
                    % Assume it is a real spike and store
                    spike_count = spike_count + 1;           
                    % now that we have the max index between -17 and +17 timesteps of
                    % the threshold we redefine put_spk_data to be centered on the
                    % max
                    max_voltage_idx = thr_cr_samp + max_timestamp_idx - 17;
                    max_put_spk_data = eeg_data(:, max_voltage_idx + (-17:17));
                    spike_mat(:,1,spike_count)     = max_voltage_idx/samp_rate*10^6; % Save spike timestamp for all channels - in microseconds - added by IV
                    spike_mat(:,2,spike_count)     = spread; % Save spread of the spike
                    spike_mat(:,3:end,spike_count) = max_put_spk_data;
%                 end

        end

        thr_cr_idx  = find( thr_cr_numind>(thr_cr_samp+35), 1, 'first' );

        if isempty( thr_cr_idx )
            break
        end

    end

end

spike_mat(:,:,spike_count+1:end) = []; % Remove unfilled matrices from spike_mat


%             % Check if the middle 5 spike samples are above the baseline
%             % value (can change to actual spike in the future)
%             spike_window = put_spk_data(max_ch, 15:19);
%             mean_spike_window = mean(spike_window); %mean voltage in the spike window 
%             is_above_baseline = abs(mean_spike_window) > abs(baseline_range);

