function [waveforms] = dat2wfs(read_dir, mapping, trial_num, trial_it, data, write_dir, data_idx, data_idx_2) 

% dat2wfs makes .dat files into arrays of voltages in the specific order
% they were recorded in and converts to waveforms using extract_wfs which
% can then be added to spatData and used in analysis. 

% INPUTS: trial directory, mapping, trial number in the folder, trial number on spatData, spatData (so 
% it can find spike times and cell IDs for a specific trial) 


% TO DO: 
% 1. I want this to make the wfs per channel also- incase i need more/want 
%    to keep channel identities  
%   tet_index is a way to keep the channel IDs and store these in the final
%   matrix - could do this post hoc 


    if nargin < 5
       write_dir = read_dir;
    end

    % load spatData and necessary data from it.   
    load(data, 'spatData');

    if ~isempty(data_idx)
        spatData = spatData (data_idx:data_idx_2, :);
    end

    %load spike times per cell in spatData 
    SpkTs = cell(size(spatData, 1), 1);
    for itC = 1:size(spatData,1)
        SpkTs{itC} = spatData.SpkTs{itC, trial_it};
    end
    %load number of spikes per cell per trial
    nSpks = spatData.nSpks(:,trial_it);
    
    
    tic;
    datFileName = dir(strcat(read_dir, '/*.dat'));
    folder = datFileName.folder;
    name = datFileName(trial_num).name
    openDat = fopen(strcat(folder,'/',name));
    
    gain = 1000; % adapting scalemax not gain 
    
    voltages = fread(openDat, [32 inf], '32*int16', 64);   % When using 'skip' format, read-in multiplier specifies N read-in format chunks, skip multiplier specifies N bytes.   
    voltages = voltages./(2^15) .* 1.5 ./gain .* -1 .* 10^6;  % TW - I rewrote this for clarity but I think it was correct already. divide by bit resolution? and multiply by voltage on scope 1.5 V - divided by gain 1000 - change to mV * 10^6 and invert trace - better for tint 
    
    % TW - Median subtract common mode noise
    voltages = voltages - median(voltages,1);
    
    % Re-order channels
    ch_reorder_ind = load(mapping).a; % variable named 'a' is totally arbitrary and can be changed in the future.
    voltages = voltages( ch_reorder_ind, : );  % TW - same function but ch order easier to read like this.
    
    %filter the data with a butterworth filter
    voltages = voltages.'; %filter works per column so it needs to be transposed 
    [b,a]    = butter(3, [300 7000]/24000, 'bandpass');
    voltages = filtfilt(b,a, voltages);

    %group voltage data into tetrodes for either single or multishank 
    if strcmp(mapping, 'map_multi_opp_plugtop_final.mat')
        % take one side vertically as a tetrode 
        tet_index = zeros(8, 4); % create an 8x4 matrix of zeros
        tet_index(1:2:end, :) = repmat(1:8:25, 4, 1).' + repmat(0:2:6, 4, 1); % fill odd-numbered rows with odd numbers
        tet_index(2:2:end, :) = repmat(2:8:26, 4, 1).' + repmat(0:2:6, 4, 1); % fill even-numbered rows with even numbers    
        %overlap tetrodes are the bottom 4 contacts which are closer to each
        %other 
        overlap = repmat (0:8:24,4,1).' + repmat (5:8,4,1);
        tet_index = [tet_index; overlap]; % + 1 overlap per octrode 
    elseif strcmp(mapping, 'map_single_opp_plugtop_final.mat') || strcmp(mapping, 'map_single_same_plugbottom_final.mat')
        tet_index = 1:32;
        tet_index = reshape(tet_index,4,8).';
        overlap_rows = repmat (0:8:24,4,1).' + repmat (3:6,4,1);
        % Add the overlapping rows to the matrix
        tet_index = [tet_index; overlap_rows];
    end 
    

    % loop through the tetrodes and make wf_mat for the current tetrode
    % if the tetrode has a cell accoring to spatData 
        
%     num_tets = size(tet_index,1); 
    waveforms = cell(height(spatData),1); 
    
    cell_tet_num = zeros(height(spatData), 1);
    for itD = 1:height(spatData) % get a list of tetrodes that have cells on them
        cell_tet_num(itD) = str2double(extractBetween(spatData.cellID(itD),'t','c'));
    end
    
    for jj = 1:height(spatData)
        for ii = cell_tet_num'
            if cell_tet_num(jj) == ii %if tetrode is on the list it runs extract_wfs 
                % Collect waveforms 
                [wf_mat] = extract_wfs( voltages(:, tet_index(ii,:) ), 48000, SpkTs, nSpks, jj); 
                % Add wfs to cell array
                waveforms{jj,1} = wf_mat;  
            end
        end
    end
    
    % save the waveforms for that trial to the write directory for that
    % experiment - if you weren't running this through collect_wfs
%     save(fullfile(write_dir, [extractBefore(name, '.') '_wfs.mat']), 'waveforms')

    % Find the index of the channel with the maximum amplitude waveform on
%     % example cell 1
%     wf_mat = exp_wfs(1,1);
%     wf_mat = wf_mat{1}; %find cell 1 trial 1 
%     [~, max_channel] = max(max(abs(wf_mat(:,:,1))));
%     % Plot all the waveforms for the channel with the maximum amplitude waveform
%     figure;
%     hold on;
%     for ii = 1:size(wf_mat,3)
%         plot(wf_mat(:,max_channel,ii));
%     end
%     xlabel('Sample Index');
%     ylabel('Amplitude');
%     title(['Waveforms for Cell 1 - Channel ' num2str(max_channel)]);
        
    toc;
end

function [wf_mat] = extract_wfs(tetrode_data,samp_rate,SpkTs, nSpks, curr_cell)
% Extract waveforms (2ms window) from all channels from corresponding
% timestamps found in spatData. Needs to have a scalemax of 200uV. 
%
% INPUT:
% tetrode_data: 4*2^19 matrix of voltage data from one file and one tetrode
% samp_rate: sample rate 
% spkTs: spike times per cell found on spatData
% nSpks: number of spikes per cell from spatData 
% curr_cell: current cell on spatData

% OUTPUT:
% wf_ mat: is a 3D array which contains for each cell all waveforms in 
% matrices of 100 samples by 4 tetrodes per spike time. 

    % Format of data should be row=ch, col=time.
    if size(tetrode_data,1)>4
        tetrode_data = tetrode_data';
    end

    % collect spkTs and add the wfs per tetrode to a matrix spike_mat

    % Preallocate waveform matrix
    wf_mat = zeros(97, 4, length(SpkTs));
    scalemax = 200; 
    % Extract waveforms for each spike time
    for ii = 1:4
        tet_spkTs = SpkTs{curr_cell, 1};
        for jj = 1:nSpks(curr_cell)
            if tet_spkTs(jj) - 0.00025 < 0
                continue;
            end
            wf_window = round(samp_rate*(tet_spkTs(jj) - 0.00025)):round(samp_rate*(tet_spkTs(jj) + 0.00175));
            wf = tetrode_data(ii, wf_window);
            %wf = wf - mean(wf); % incase i need to center the waveforms
            %because of a DC offset 
            wf = (wf ./ scalemax)*128; % not sure if this is good 
            wf_mat(:, ii, jj) = wf.';
        end
    end
end





  

    
