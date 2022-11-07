 function [spike_mat,spike_count] = extract_spikes_TW(tetrode_data,samp_rate,spikeThreshold)
% Extract spike data (1ms window) from all channels and extract corresponding
% timestamp. Spikes are considered voltages which cross the DACQ threshold
% (70.3 microvolts).
%
% INPUT:
% tetrode_data: 4*2^19 matrix of voltage data from one file and one tetrode
%
% OUTPUT:
% spike_count: total number of spikes recorded
% spike_mat: matrix containing spike time (row 1) and corresponding
%   voltages across all channels.

% Create empty matrix to ensure correct timestamp stored

% % If its not the first file with Bonsai delay then will be a full file 
% if nargin < 4 
%     empty_matrix = nan(length(tetrode_data),32)'; % size of one file including 16 channels % Iv changed 16 to 32 and 2^19 to samples (length of voltages matrix)  
% else
%     empty_matrix = nan(size(tetrode_data,1),32)'; % size of fill with delay removed
% end

% Format of data should be row=ch, col=time.
if size(tetrode_data,1)>4
    tetrode_data = tetrode_data';
end

% Insert tetrode data at correct position in the empty matrix % consider
% getting rid of iterations - since this is a tara code thing - will remove
% the need to loop over this function in dat2spikes and produce ony one
% spike_mat 
% if iteration == 1
%     empty_matrix(1:4,:) = tetrode_data;
% elseif iteration == 2
%     empty_matrix(5:8,:) = tetrode_data;
% elseif iteration == 3
%     empty_matrix(9:12,:) = tetrode_data;
% elseif iteration == 4
%     empty_matrix(13:16,:) = tetrode_data;
% elseif iteration == 5
%     empty_matrix(17:20,:) = tetrode_data;
% elseif iteration == 6
%     empty_matrix(21:24,:) = tetrode_data;
% elseif iteration == 7
%     empty_matrix(25:28,:) = tetrode_data;
% elseif iteration == 8
%     empty_matrix(29:32,:) = tetrode_data;
% end


thr_cr         = max(tetrode_data,[],1) > spikeThreshold;  % 'thr_cr' is a logical index of where threshold is crossed
thr_cr_numind  = find( thr_cr );                  % 'thr_cr_numind' is a numerical index for where threshold crossed

possible_count    = sum(thr_cr);              % Maximum possible no. spikes. TW - not great like this, prob x10 bigger than actual, but actual N spikes difficult to predict. Should probably use total contiguous crossings.
spike_mat         = nan(4,51,possible_count); % Matrix to store spike data
nSamp             = size( tetrode_data, 2 );

upperThreshold    = 250; % ignore 'spikes' bigger than 500uv % IV changed to 200

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
        [max_val,max_ind] = max(put_spk_data,[],'all');
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



















