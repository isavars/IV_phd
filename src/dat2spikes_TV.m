% dat2spikes makes .dat files into arrays of voltages in the specific order they were recorded in and converts to spikes using tara's extract spikes funciton and saves them as a
% .mat file to be used by Tara's tint conversion pipeline to be outputed as
% tint files. 

% The input to this function is a trial directory, mapping
function [] = dat2spikes_TV(read_dir, mapping, write_dir)
    if nargin < 3
       write_dir = read_dir;
    end

    tic;
    datFileName = dir(strcat(read_dir, '/*.dat'));
    folder = datFileName.folder;
    name = datFileName.name;
    openDat = fopen(strcat(folder,'/',name));
    
    
    voltages = fread(openDat, [32 inf], '32*int16', 64);   % When using 'skip' format, read-in multiplier specifies N read-in format chunks, skip multiplier specifies N bytes.   
    voltages = voltages./(2^15) .* 1.5 ./1000 .* -1 .* 10^6;                  % TW - I rewrote this for clarity but I think it was correct already. divide by bit resolution? and multiply by voltage on scope 1.5 V - divided by gain 1000 - change to mV * 10^6 and invert trace - better for tint 
    
    % TW - Median subtract common mode noise
    voltages = voltages - median(voltages,1);
    
    % Re-order channels
    ch_reorder_ind = mapping;
    voltages = voltages( ch_reorder_ind, : );  % TW - same function but ch order easier to read like this.
    
    %filter the data with a butterworth filter
    voltages = voltages.'; %filter works per column so it needs to be transposed 
    [b,a]    = butter(3, [300 7000]/24000, 'bandpass');
    voltages = filtfilt(b,a, voltages);
    %     voltages = voltages([3000:600.0625*48000],:);
    
    %get data into exact shape for tara's code 
    %group voltage data into tetrodes 
    tet_index = (1:4) + (0:4:28)';    
    %overlap tetrodes     
    
    % run voltage traces through extract_spikes
    
    duration = length(voltages)/48000 ; %trial duration for timestamps divide #samples by sampling rate     
    % loop through the tetrodes and make spike_mat and spike_count for each to
    % be used by makeTetrodes
        
    num_tets = size(tet_index,1); % chance to 12 if you do the overlap
    name = datFileName.name;
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
    %     tetrode9 = cell2mat(new_tets(1,:)');
    %     tetrode10 = cell2mat(new_tets(2,:)');
    %     tetrode11 = cell2mat(new_tets(3,:)');
    %     tetrode12 = cell2mat(new_tets(4,:)');    
    
    % Sum all spikes in trial
    num_spikes = sum(cell2mat(num_spikes),2);
    
    %cut trtial to time 
    
    trial_duration = round(duration);
    
    [tetrode1,num_spikes(1)] =  cut_trial_to_time(tetrode1,trial_duration,num_spikes(1));
    [tetrode2,num_spikes(2)] =  cut_trial_to_time(tetrode2,trial_duration,num_spikes(2));
    [tetrode3,num_spikes(3)] =  cut_trial_to_time(tetrode3,trial_duration,num_spikes(3));
    [tetrode4,num_spikes(4)] =  cut_trial_to_time(tetrode4,trial_duration,num_spikes(4));
    [tetrode5,num_spikes(5)] =  cut_trial_to_time(tetrode5,trial_duration,num_spikes(5));
    [tetrode6,num_spikes(6)] =  cut_trial_to_time(tetrode6,trial_duration,num_spikes(6));
    [tetrode7,num_spikes(7)] =  cut_trial_to_time(tetrode7,trial_duration,num_spikes(7));
    [tetrode8,num_spikes(8)] =  cut_trial_to_time(tetrode8,trial_duration,num_spikes(8));
    
    name = datfileName.name;
    name = extractBefore(name, ".");
    trialInfo.finalTrialName = name;
    trialInfo.writeDir = write_dir;
    
    setFileName = dir(strcat(read_dir, '/*.set'));
    fopen([setFileName.name '.set']);
    % Read file %
    [key value] = textread([setFileName.name], '%s %[^\n]');
    txt = [cat(1,key) cat(1,value)];
    trialInfo.trial_time = getValue(txt, 'trial_time');
    trialInfo.trial_date = getValue(txt, 'trial_date');
    trialInfo.duration = getValue(txt, 'duration');
        
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
    
    toc;
end
%     save spikes and spike count from extract spikes output/ other variabes you want to save
   




    