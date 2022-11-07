% dat2spikes makes .dat files into arrays of voltages in the specific order they were recorded in and converts to spikes using tara's extract spikes funciton and saves them as a
% .mat file to be used by Tara's tint conversion pipeline to be outputed as
% tint files. 


    tic;

    fileName = dir('**/*.dat');
    name = fileName.name;
    openDat = fopen(name);

%     voltages = fread(openDat, [64 48000*30], 'int16');   % TW - just use one variable for memory efficiency
%     voltages = reshape(voltages, 64, []);
%     voltages = voltages(1:32,:);                       %CM32 take out the empty bits 
    
    voltages = fread(openDat, [32 inf], '32*int16', 64);   % When using 'skip' format, read-in multiplier specifies N read-in format chunks, skip multiplier specifies N bytes.   
    voltages = voltages./(2^15) .* 1.5 ./1000 .* -1 .* 10^6;                  % TW - I rewrote this for clarity but I think it was correct already. divide by bit resolution? and multiply by voltage on scope 1.5 V - divided by gain 1000 - change to mV * 10^6 and invert trace - better for tint 
    
    % TW - Median subtract common mode noise
    voltages = voltages - median(voltages,1);

    %reorganize row order to fit map for each type of probe 
%     % Re-ordering based on CM32 Buzsaki 4x8 plugged in with NN and Omnetics on opposite sides 
%     ch_reorder_ind = [6 2 5 29     27 3 4 28 ...
%                       30 8 1 7     31 25 32 26 ...
%                       9 19 10 16   24 18 23 17 ...
%                       15 11 20 12  14 22 21 13];
%     BUZ OPP plugTop
%     ch_reorder_ind = [11 15 12 20   22 14 13 21 ...
%                       19 9 16 10   18 24 17 22 ...
%                       8 30 7 1      25 31 26 32  ...
%                       2 6 29 5      3 27 28 4];
    % BUZ SAME plugTop
%     ch_reorder_ind = [27 31 28 4   6 30 29 5 ...
%                       3 25 32 26   2 8 1 6 ...
%                       24 14 23 17      9 15 10 16  ...
%                       18 22 13 21      19 11 12 20];
    %BUZ SAME
%     ch_reorder_ind = [22 18 21 13   11 19 20 12 ...
%                   14 24 17 23   15 9 16 10 ...
%                   25 3 26 32      8 2 7 1  ...
%                   31 27 4 28      30 6 5 29];

%     %Tom's Map
%         ch_reorder_ind = [1 2 6 7   8 9 15 16 ...
%                       3 4 5 10   11 12 13 14 ...
%                       19 20 21 22  23 28 29 30 ...  
%                       17 18 24 25  26 27 31 32];
%     %Tom's Map - reshuffled 
%         ch_reorder_ind = [ 9 15 16 6  2 1 8 7 ...
%                           4 12 5 11   3 13 10 14 ...
%                           21 29 30 22   28 20 19 23 ...
%                           17 27 18 24   26 32 31 25];

    % Re-ordering based on CM32 1x32 OPP plugTop
% 
    ch_reorder_ind = [8 9 16 1   7 19 30 10 ...
                  15 2 25 20   29 24 14 3 ...
                  26 21 28 23      13 4 32 22  ...
                  27 17 12 5      31 11 6 18];

   % Re-ordering based on CM32 1x32 SAME plugTop

    ch_reorder_ind = [24 25 32 17   23 3 14 26 ...
                    31 18 9 4   13 8 30 19 ...
                    10 5 12 7   29 20 16 6 ...
                    11 1 28 21  15 27 22 2];

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
%     tet_index = (1:4) + (2:8:30)';
    
    
    % run voltage traces through extract_spikes
    
     duration = length(voltages)/48000 ; %trial duration for timestamps divide #samples by sampling rate 
     duration_us = duration*10^6;
%     sample = 1/48000; % one sample in seconds 
%     sample_uV = 1/48000*10^6; %length of one sample in microseconds 
%     timestamps = 0:sample_uV:duration_us-1;
    
    % loop through the tetrodes and make spike_mat and spike_count for each to
    % be used by makeTetrodes
        
        num_tets = size(tet_index,1); % chance to 12 if you do the overlap
        name = fileName.name;
        fileNames = {extractBefore(name, ".")};
        new_tets = cell(num_tets,numel(fileNames)); %ive changed these to 8 because I think they correspond to number of tetrodes. 
        num_spikes = cell(num_tets,numel(fileNames));
        num_spikes(:,:) = {0};
        jj =1; 
    
%         while jj <= numel(fileNames)
            for ii = 1:num_tets
                % Collect spike data
                [spike_mat,spike_count] = extract_spikes_TW( voltages(:, tet_index(ii,:) ), 48000, 60); 
                % Increment spike count on channel
                num_spikes{ii,jj} = spike_count;   
                % Reshape matrix 
                [final_mat] = reshape_spike_mat(spike_mat);
                % Add spike data to cell array
                new_tets{ii,jj} = final_mat;  
            end
    %         jj = jj +1;
%         end
    
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
 
    toc;
    %     save spikes and spike count from extract spikes output/ other variabes you want to save
       
    
    


    