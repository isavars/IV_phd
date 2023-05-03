function [spike_mat,spike_count] = extract_wfs(tetrode_data,samp_rate,spkTs)
% Extract waveforms (2ms window) from all channels from corresponding
% timestamps found in spatData. 
%
% INPUT:
% tetrode_data: 4*2^19 matrix of voltage data from one file and one tetrode
% samp_rate: sample rate 
% spkTs: spike times from cut files found on spatData 

% OUTPUT:
% waveforms and waveform means: in the same format that spatData is organized 


% Format of data should be row=ch, col=time.
if size(tetrode_data,1)>4
    tetrode_data = tetrode_data';
end
% total_nSpks - get from spatData
all_spikes    = sum(total_nSpks);              % get this from spatData
spike_mat         = nan(4,51,all_spikes); % Matrix to store spike data -  is this where I can just take 100 samples instead of 50? 
nSamp             = size( tetrode_data, 2 );

spike_count = 0;
thr_cr_idx  = find( thr_cr_numind>10, 1, 'first' );    % 'thr_cr_idx' is the iterator over ''thr_cr_numind' for the loop below. Get the earliest >sample 10, so as to allow pre-crossing samples.

% Loop through crossings and extract spikes
while thr_cr_idx <= length(thr_cr_numind)
    
    thr_cr_samp = thr_cr_numind(  thr_cr_idx  );  % 'thr_cr_samp' is the current number of the sample where threshold is crossed.
    
    if thr_cr_samp > nSamp-40  % If we have reached the end of the data, stop collecting spikes
        break
    else 

        put_spk_data = tetrode_data( :, thr_cr_samp + (-10:39) );  % 'put_spk_data' = putative spike data, putative as hasn't undergone artifact filtering yet.
        
        % Get the necessary info for testing artifact rejection
        [max_val,max_ind] = max(put_spk_data,[],'all', 'linear');
        [max_ch,~]        = ind2sub( size(put_spk_data), max_ind );
        first_samp        = put_spk_data( max_ch, 1);
        late_samps        = put_spk_data( max_ch, 41:50);

        if (first_samp < max_val*0.78) &&  (first_samp > -max_val*0.78)  &&  ~any(late_samps>max_val*0.33)  &&  max_val<upperThreshold % This is artefact rejection test
            
            % Assume is a real spike and store
            spike_count                    = spike_count + 1;
            spike_mat(:,1,spike_count)     = thr_cr_samp/samp_rate*10^6; % Save spike timestamp for all channels - in microseconds - added by IV
            spike_mat(:,2:end,spike_count) = put_spk_data; %wont add the voltages to the matrix without the interpolation - its not actually needed but it might not cause any harm.                          

        end

        thr_cr_idx  = find( thr_cr_numind>(thr_cr_samp+50), 1, 'first' );

        if isempty( thr_cr_idx )
            break
        end

    end

end

spike_mat(:,:,spike_count+1:end) = []; % Remove unfilled matrices from spike_mat



















