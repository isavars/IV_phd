% dat2spikes makes .dat files into arrays of voltages in the specific order they were recorded in and converts to spikes using tara's extract spikes funciton and saves them as a
% .mat file to be used by Tara's tint conversion pipeline to be outputed as
% tint files. 

% The input to this function is a trial directory, mapping
function [] = dat2spikes_IV_old(read_dir, mapping, trial_num, write_dir)
    if nargin < 4
       write_dir = read_dir;
    end

    tic;
    datFileName = dir(strcat(read_dir, '/*.dat'));
    folder = datFileName.folder;
    name = datFileName(trial_num).name;
    openDat = fopen(strcat(folder,'/',name));
    
    
    voltages = fread(openDat, [32 inf], '32*int16', 64);   % When using 'skip' format, read-in multiplier specifies N read-in format chunks, skip multiplier specifies N bytes.   
    voltages = voltages./(2^15) .* 1.5 ./1000 .* -1 .* 10^6;                  % TW - I rewrote this for clarity but I think it was correct already. divide by bit resolution? and multiply by voltage on scope 1.5 V - divided by gain 1000 - change to mV * 10^6 and invert trace - better for tint 
    
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
    
%     duration = length(voltages)/48000 ; %trial duration for timestamps divide #samples by sampling rate    

    % loop through the tetrodes and make spike_mat and spike_count for each to
    % be used by makeTetrodes
        
    num_tets = size(tet_index,1); % change to 12 if you do the overlap
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
        
    % Concatenate info from each file - this is what tara needs to do but I
    % sample continously so maybe not needed? is she combining each tetrode for
    % a full experoment or is this a per file thing? 
    
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
    
    % Sum all spikes in trial
    num_spikes = sum(cell2mat(num_spikes),2);
    
    
    %read in set file 
    name = extractBefore(name, ".");
    trialInfo.finalTrialName = name;
    display(write_dir)
    trialInfo.writeDir = write_dir;
    
    setFileName = dir(strcat(read_dir, '/*.set'));
    folder = setFileName.folder;
    name = setFileName(trial_num).name;
    filepath = strcat(folder, '/', name);
    fopen([filepath '.set']);
    % Read file %
    [key value] = textread([filepath], '%s %[^\n]');
    txt = [cat(1,key) cat(1,value)];
    trialInfo.trial_time = getValue(txt, 'trial_time');
    trialInfo.trial_date = getValue(txt, 'trial_date');
    trialInfo.duration = getValue(txt, 'duration');
    
    %cut trtial to time 
    trial_duration = str2double(trialInfo.duration);
    [tetrode1,num_spikes(1)] =  cut_trial_to_time(tetrode1,trial_duration,num_spikes(1));
    [tetrode2,num_spikes(2)] =  cut_trial_to_time(tetrode2,trial_duration,num_spikes(2));
    [tetrode3,num_spikes(3)] =  cut_trial_to_time(tetrode3,trial_duration,num_spikes(3));
    [tetrode4,num_spikes(4)] =  cut_trial_to_time(tetrode4,trial_duration,num_spikes(4));
    [tetrode5,num_spikes(5)] =  cut_trial_to_time(tetrode5,trial_duration,num_spikes(5));
    [tetrode6,num_spikes(6)] =  cut_trial_to_time(tetrode6,trial_duration,num_spikes(6));
    [tetrode7,num_spikes(7)] =  cut_trial_to_time(tetrode7,trial_duration,num_spikes(7));
    [tetrode8,num_spikes(8)] =  cut_trial_to_time(tetrode8,trial_duration,num_spikes(8));
    [tetrode9,num_spikes(9)] =  cut_trial_to_time(tetrode9,trial_duration,num_spikes(9));
    [tetrode10,num_spikes(10)] =  cut_trial_to_time(tetrode10,trial_duration,num_spikes(10));
    [tetrode11,num_spikes(11)] =  cut_trial_to_time(tetrode11,trial_duration,num_spikes(11));
    [tetrode12,num_spikes(12)] =  cut_trial_to_time(tetrode12,trial_duration,num_spikes(12));
    %% Generate Tint .tet header
    [genericTetHeader] = TintTetrode_header(trialInfo);%replaced trialInfo,DACQinfo with duration
    [tet1FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,1);
    [tet2FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,2);
    [tet3FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,3);
    [tet4FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,4);
    [tet5FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,5);
    [tet6FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,6);
    [tet7FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,7);
    [tet8FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,8);
    [tet9FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,9);
    [tet10FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,10);
    [tet11FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,11);
    [tet12FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,12);
    %% Append header and create tetrode file
    fprintf('Generating tetrode files....\n')
    make_tint_file(tetrode1,tet1FileName,1,trialInfo); % removed ,DACQinfo
    make_tint_file(tetrode2,tet2FileName,2,trialInfo);
    make_tint_file(tetrode3,tet3FileName,3,trialInfo);
    make_tint_file(tetrode4,tet4FileName,4,trialInfo);
    make_tint_file(tetrode5,tet5FileName,5,trialInfo);
    make_tint_file(tetrode6,tet6FileName,6,trialInfo);
    make_tint_file(tetrode7,tet7FileName,7,trialInfo);
    make_tint_file(tetrode8,tet8FileName,8,trialInfo);
    make_tint_file(tetrode9,tet9FileName,9,trialInfo);
    make_tint_file(tetrode10,tet10FileName,10,trialInfo);
    make_tint_file(tetrode11,tet11FileName,11,trialInfo);
    make_tint_file(tetrode12,tet12FileName,12,trialInfo);
    toc;
end
%     save spikes and spike count from extract spikes output/ other variabes you want to save
   




    