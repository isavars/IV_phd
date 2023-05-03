% dat2wfs makes .dat files into arrays of voltages in the specific order
% they were recorded in and converts to waveforms using extract_wfs which
% can then be added to spatData and used in analysis. 

% The input to this function is a trial directory, mapping - want this to
% loop through spatData and find all the relevant .dat files using rat numbers and
% dates in the folders in the big computer - ask TV 

function [] = dat2wfs(read_dir, mapping, trial_num, write_dir)
    
    % load spatData and necessary data from it. 

    load("spatData_probes.mat",'spatData');

    if nargin < 4
       write_dir = read_dir;
    end
    
    tic;
    datFileName = dir(strcat(read_dir, '/*.dat'));
    folder = datFileName.folder;
    name = datFileName(trial_num).name;
    openDat = fopen(strcat(folder,'/',name));
    
    %gain adjustment per file - method for getting rat id can switch to
    %getting it from spatData

    load('gains_thesis_data.mat', 'gains') % trying to avoid clipping completely here so do I just set new very conservative gains? - ask TW 
    curr_rat = extractBetween(read_dir,'r','/');
    curr_rat = str2double(curr_rat{2});
    curr_rat_idx = find(gains.rat_ID == curr_rat);
    gain = gains.gain(curr_rat_idx); 
    
    voltages = fread(openDat, [32 inf], '32*int16', 64);   % When using 'skip' format, read-in multiplier specifies N read-in format chunks, skip multiplier specifies N bytes.   
    voltages = voltages./(2^15) .* 1.5 ./gain .* -1 .* 10^6;                  % TW - I rewrote this for clarity but I think it was correct already. divide by bit resolution? and multiply by voltage on scope 1.5 V - divided by gain 1000 - change to mV * 10^6 and invert trace - better for tint 
    
    % TW - Median subtract common mode noise
    voltages = voltages - median(voltages,1);
    
    % Re-order channels
    ch_reorder_ind = load(mapping).a; % variable named 'a' is totally arbitrary and can be changed in the future.
    voltages = voltages( ch_reorder_ind, : );  % TW - same function but ch order easier to read like this.
    
    %filter the data with a butterworth filter
    voltages = voltages.'; %filter works per column so it needs to be transposed 
    [b,a]    = butter(3, [300 7000]/24000, 'bandpass');
    voltages = filtfilt(b,a, voltages);
    %     voltages = voltages([3000:600.0625*48000],:);
    
    %get data into exact shape for tara's code 
    %group voltage data into tetrodes 
    % take one side vertically as a tetrode 
    tet_index = zeros(8, 4); % create an 8x4 matrix of zeros
    tet_index(1:2:end, :) = repmat(1:8:25, 4, 1).' + repmat(0:2:6, 4, 1); % fill odd-numbered rows with odd numbers
    tet_index(2:2:end, :) = repmat(2:8:26, 4, 1).' + repmat(0:2:6, 4, 1); % fill even-numbered rows with even numbers

    %overlap tetrodes are the bottom 4 contacts which are closer to each
    %other 
    overlap = repmat (0:8:24,4,1).' + repmat (5:8,4,1);
    tet_index = [tet_index; overlap]; % + 1 overlap per octrode 
    
    % run voltage traces through extract_spikes

    % loop through the tetrodes and make spike_mat and spike_count for each to
    % be used by makeTetrodes
        
    num_tets = size(tet_index,1); 
%     name = datFileName(trial_num).name;
    fileNames = {extractBefore(name, ".")};
    new_tets = cell(num_tets,numel(fileNames)); %ive changed these to 8 because I think they correspond to number of tetrodes. 
    num_spikes = cell(num_tets,numel(fileNames));
    num_spikes(:,:) = {0};
    
    for ii = 1:num_tets
        % Collect spike data
        [spike_mat,spike_count] = extract_spikes_TW( voltages(:, tet_index(ii,:) ), 48000, 60); 
        % Increment spike count on channel
        num_spikes{ii,1} = spike_count;   
        % Reshape matrix 
        [final_mat] = reshape_spike_mat(spike_mat);
        % Add spike data to cell array
        new_tets{ii,1} = final_mat;  
    end

    % Concatenate info from each file - do i need to make these? maybe I
    % just want the wfs per channel. 
    
    tetrode1 = cell2mat(new_tets(1,:)');
    tetrode2 = cell2mat(new_tets(2,:)');
    tetrode3 = cell2mat(new_tets(3,:)');
    tetrode4 = cell2mat(new_tets(4,:)');
    tetrode5 = cell2mat(new_tets(5,:)');
    tetrode6 = cell2mat(new_tets(6,:)');
    tetrode7 = cell2mat(new_tets(7,:)');
    tetrode8 = cell2mat(new_tets(8,:)');
    tetrode9 = cell2mat(new_tets(9,:)');
    tetrode10 = cell2mat(new_tets(10,:)');
    tetrode11 = cell2mat(new_tets(11,:)');
    tetrode12 = cell2mat(new_tets(12,:)');   

    % might need to add the cut to time bit here - sleep trials are going
    % to cause a real problem because the cut files are joined and they are
    % trimmed down to sws.- might have to work from the scan files instead?
        
    toc;
end
%     make waveform means and save waveforms and waveform means per channel
%     so they can be accessed when making get spat data or appended easilly
%     to it. - alternatively they never get added and its only used when
%     running the waveform PCA. The channel map needs to be somewhere
%     though. 
