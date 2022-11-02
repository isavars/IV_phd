 function [spike_mat,spike_count] = extract_spikes(tetrode_data,iteration,timestamps,spikeThreshold)
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

%% Create empty matrix to ensure correct timestamp stored

upperThreshold = 500; % ignore 'spikes' bigger than 500uv % IV changed to 200

% % If its not the first file with Bonsai delay then will be a full file 
% if nargin < 4 
    empty_matrix = nan(length(tetrode_data),32)'; % size of one file including 16 channels % Iv changed 16 to 32 and 2^19 to samples (length of voltages matrix)  
% else
%     empty_matrix = nan(size(tetrode_data,1),32)'; % size of fill with delay removed
% end
tetrode_data = tetrode_data';

% Insert tetrode data at correct position in the empty matrix % consider
% getting rid of iterations - since this is a tara code thing - will remove
% the need to loop over this function in dat2spikes and produce ony one
% spike_mat 
if iteration == 1
    empty_matrix(1:4,:) = tetrode_data;
elseif iteration == 2
    empty_matrix(5:8,:) = tetrode_data;
elseif iteration == 3
    empty_matrix(9:12,:) = tetrode_data;
elseif iteration == 4
    empty_matrix(13:16,:) = tetrode_data;
elseif iteration == 5
    empty_matrix(17:20,:) = tetrode_data;
elseif iteration == 6
    empty_matrix(21:24,:) = tetrode_data;
elseif iteration == 7
    empty_matrix(25:28,:) = tetrode_data;
elseif iteration == 8
    empty_matrix(29:32,:) = tetrode_data;
end

%% Store spike data and associated timestamp across all channels

cross_thresh = tetrode_data > spikeThreshold;
possible_count = sum(sum(cross_thresh)); % Maximum possible no. spikes 

spike_mat = nan(4,51,possible_count); % Matrix to store spike data

spike_count = 0;
idx = 1;
spike_val = nan(1,possible_count);

% Get 9 samples before spike + spike/threshold crossing + 40 samples after
while idx <= length(empty_matrix) %changing numel to length
    [n,m] = ind2sub(size(empty_matrix),idx);
    n = n - 4*(iteration-1);
    if m < 10
        idx = idx + minus(10,m);
    elseif m >= 10 && empty_matrix(idx) > spikeThreshold && empty_matrix(idx) < upperThreshold % ignore 'spikes' bigger than upperThreshold uV
        spike_count = spike_count + 1; % Increment spike count 
        spike_val(spike_count) = empty_matrix(idx); % Get voltage value          
        if (m - 6) <= 0 || (m + 26) > size(tetrode_data,2)
            spike_count = spike_count - 1;
            if spike_count > 0;  spike_val(spike_count) = nan; end 
            if(m + 26) > size(tetrode_data,2); break; end
        else
            % Check if spike is dodgy/noise (Baseline = 0)
            % Based on artefact rejection criteria from DACQ manual (page 44-45)
            first_sample = tetrode_data(n,m-6);
            last_20_samples = tetrode_data(n,m+6:1:m+25);
            % Reject if first sample well above or below baseline (>78% or <-78% spike)
            if first_sample > (spike_val(spike_count)*0.78) || first_sample < -1*(spike_val(spike_count)*0.78)
                 spike_count = spike_count - 1;
            % Also reject if last 20 samples are > 33% above baseline
            elseif sum(last_20_samples > (0.33*spike_val(spike_count))) > 1
                 spike_count = spike_count - 1;
            else % Assume is a real spike
                spike_mat(:,1,spike_count) = repelem(timestamps(idx),4); % Save spike timestamp for all channels
                % Interpolate
                interpolateSpike = tetrode_data(:,m-6:1:m+25);
                samplePoints = 1:32;
                queryPoints = linspace(1,32,50);
                spike_mat(:,2:end,spike_count)= interp1(samplePoints,interpolateSpike',queryPoints,'linear')'; %wont add the voltages to the matrix without the interpolation - its not actually needed but it might not cause any harm.                          
            end
            idx = sub2ind(size(empty_matrix),1,m+26); % Jump to end of 1ms window
        end                   
    else
          idx = idx + 1;
    end

end

spike_mat(:,:,spike_count+1:end) = []; % Remove unfilled matrices from spike_mat



















