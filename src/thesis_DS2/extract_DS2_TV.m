function [spike_mat,spike_count] = extract_DS2_TV(eeg_data,samp_rate,spikeThreshold)
% Extract spike data (50ms window) from all channels and extract corresponding
% timestamp. Spikes are considered voltages which cross the DACQ threshold
% (1.14 V).
% 140ms - 35 samples
% INPUT:
% eeg_data: 32*samples matrix of voltage data from all channels in sleep
% trial
%
% OUTPUT:
% spike_count: total number of dentate spikes recorded
% spike_mat: matrix containing spike time (row 1) and corresponding
%   voltages across all channels.

eeg_data_abs= abs(eeg_data);%get the absiolute value so positive and negative threshold values are considered. 

thr_cr         = max(eeg_data_abs,[],1) > spikeThreshold;  % 'thr_cr' is a logical index of where threshold is crossed
thr_cr_numind  = find( thr_cr );                  % 'thr_cr_numind' is a numerical index for where threshold crossed

possible_count    = sum(thr_cr);              % Maximum possible no. spikes. 
spike_mat         = nan(32,37,possible_count); % Matrix to store spike data
nSamp             = size( eeg_data, 2 );

upperThreshold    = 250; % ignore 'spikes' bigger than 500uv % IV changed to 200

spike_count = 0;
thr_cr_idx  = find( thr_cr_numind>15, 1, 'first' );    % 'thr_cr_idx' is the iterator over ''thr_cr_numind' for the loop below. Get the earliest >sample 10, so as to allow pre-crossing samples.

% Loop through crossings and extract spikes
while thr_cr_idx <= length(thr_cr_numind)
    
    thr_cr_samp = thr_cr_numind(  thr_cr_idx  );  % 'thr_cr_samp' is the current number of the sample where threshold is crossed.
    
    if thr_cr_samp > nSamp-20  % If we have reached the end of the data, stop collecting spikes
        break
    else 

        put_spk_data = eeg_data_abs( :, thr_cr_samp + (-17:17) );  % 'put_spk_data' = putative spike data, putative as hasn't undergone artifact filtering yet.
        raw_put_spk_data = eeg_data( :, thr_cr_samp + (-17:17) );
        % Get the necessary info for testing artifact rejection
        [max_val,max_ind] = max(put_spk_data,[],'all', 'linear');
        [max_ch, max_timestamp_idx]= ind2sub( size(put_spk_data), max_ind);
        
        first_samp        = put_spk_data( max_ch, 1);
        late_samps        = put_spk_data( max_ch, 25:35);%assuming last 5 samples should be flat
        spread_array = raw_put_spk_data(:, max_timestamp_idx);
        spread = max(spread_array) - min(spread_array);
        % needs both big spread at spike
        if (first_samp < max_val*0.78) &&  (first_samp > -max_val*0.78)  &&  ~any(late_samps>max_val*0.33)  &&  max_val<upperThreshold % This is artefact rejection test
            
            % Assume is a real spike and store
            spike_count                    = spike_count + 1;
            % now that we have the max index between -17 and +17 timesteps of
            % the threshold we redefine put_spk_data to be centered on the
            % max
            max_voltage_idx = thr_cr_samp + max_timestamp_idx - 17;
            max_put_spk_data = eeg_data(:, max_voltage_idx + (-17:17));
            spike_mat(:,1,spike_count)     = max_voltage_idx/samp_rate*10^6; % Save spike timestamp for all channels - in microseconds - added by IV
            spike_mat(:,2,spike_count)     = spread; % Save spread of the spike
            spike_mat(:,3:end,spike_count) = max_put_spk_data;
        end

        thr_cr_idx  = find( thr_cr_numind>(thr_cr_samp+35), 1, 'first' );

        if isempty( thr_cr_idx )
            break
        end

    end

end

spike_mat(:,:,spike_count+1:end) = []; % Remove unfilled matrices from spike_mat



