% dat2spikes makes .dat files into arrays of voltages in the specific order they were recorded in and converts to spikes using tara's extract spikes funciton and saves them as a
% .mat file to be used by Tara's tint conversion pipeline to be outputed as
% tint files. 


    tic;

    fileName = dir('**/*.dat');
    openDat = fopen(fileName.name);
    readDat = fread(openDat, 'int16');
    voltages = reshape(readDat, 64, []);
    voltages = voltages([1:32],:);%CM32 take out the empty bits 
    voltages = voltages./(2^15).*1.5./1000.*10^6.* -1; %divide by bit resolution? and multiply by voltage on scope 1.5 V - divided by gain 1000 - change to mV * 10^6 and invert trace - better for tint 
    clear readDat
    %reorganize row order to fit map for each type of probe 

    % Re-ordering based on CM32 Buzsaki 4x8 plugged in with NN and Omnetics on opposite sides 
    %voltages = [voltages(6,:);voltages(2,:);voltages(5,:);voltages(29,:);voltages(27,:);voltages(3,:);voltages(4,:);voltages(28,:);voltages(30,:);voltages(8,:);voltages(1,:);voltages(7,:);voltages(31,:);voltages(25,:);voltages(32,:);voltages(26,:);voltages(9,:);voltages(19,:);voltages(10,:);voltages(16,:);voltages(24,:);voltages(18,:);voltages(23,:);voltages(17,:);voltages(15,:);voltages(11,:);voltages(20,:);voltages(12,:);voltages(14,:);voltages(22,:);voltages(21,:);voltages(13,:)];
    %opposite pluging option - incase the map is wong 
    voltages = [voltages(22,:);voltages(18,:);voltages(21,:);voltages(13,:);voltages(11,:);voltages(19,:);voltages(20,:);voltages(12,:);voltages(14,:);voltages(24,:);voltages(17,:);voltages(23,:);voltages(15,:);voltages(9,:);voltages(16,:);voltages(10,:);voltages(25,:);voltages(3,:);voltages(26,:);voltages(32,:);voltages(8,:);voltages(2,:);voltages(7,:);voltages(1,:);voltages(31,:);voltages(27,:);voltages(4,:);voltages(28,:);voltages(30,:);voltages(6,:);voltages(5,:);voltages(29,:)];
    
    % Re-ordering based on CM32 1x32 plugged in with NN and Omnetics on opposite sides
    %voltages = [voltages(8,:);voltages(9,:);voltages(16,:);voltages(1,:);voltages(7,:);voltages(19,:);voltages(30,:);voltages(10,:);voltages(15,:);voltages(2,:);voltages(25,:);voltages(20,:);voltages(29,:);voltages(24,:);voltages(14,:);voltages(3,:);voltages(26,:);voltages(21,:);voltages(28,:);voltages(23,:);voltages(13,:);voltages(4,:);voltages(32,:);voltages(22,:);voltages(27,:);voltages(17,:);voltages(12,:);voltages(5,:);voltages(31,:);voltages(11,:);voltages(6,:);voltages(18,:)];
    
    %filter the data with a butterworth filter 
    [b,a] = butter(3, [300 7000]/24000, 'bandpass');
    voltages = voltages.'; %filter works per column so it needs to be transposed 
    voltages = filtfilt(b,a, voltages);
    
    %get data into exact shape for tara's code 
    
    %group voltage data into tetrodes 
    
     voltages_1 = (voltages(:,[1:4]));
     voltages_2 = (voltages(:,[5:8]));
     voltages_3 = (voltages(:,[9:12]));
     voltages_4 = (voltages(:,[13:16]));
     voltages_5 = (voltages(:,[17:20]));
     voltages_6 = (voltages(:,[21:24]));
     voltages_7 = (voltages(:,[25:28]));
     voltages_8 = (voltages(:,[29:32]));
    
     % overlap tetrodes - just in case 
    %  voltages_9 = (voltages(:,[3:6]));
    %  voltages_10 = (voltages(:,[11:14]));
    %  voltages_11 = (voltages(:,[19:22]));
    %  voltages_12 = (voltages(:,[27:30]));
    
    %create a cell array that is made up of voltages grouped by tetrodes 
    
    all_tets= {voltages_1,voltages_2,voltages_3,voltages_4,voltages_5,voltages_6,voltages_7,voltages_8}; %voltages_9,voltages_10,voltages_11,voltages_12};
    
    
    % run voltage traces through extract_spikes
    
    duration = length(voltages)/48000 ; %trial duration for timestamps divide #samples by sampling rate 
    duration_us = duration*10^6;
    sample = 1/48000; %one sample in seconds 
    sample_uV = 1/48000*10^6; %length of one sample in microseconds 
    timestamps = 0:sample_uV:duration_us-1;
    
    % loop through the tetrodes and make spike_mat andd spike_count for each to
    % be used by makeTetrodes
        
        tetNum = 8; %chance to 12 if you do the overlap
        name = fileName.name;
        fileNames = {extractBefore(name, ".")};
        new_tets = cell(8,numel(fileNames)); %ive changed these to 8 because I think they correspond to number of tetrodes. 
        num_spikes = cell(8,numel(fileNames));
        num_spikes(:,:) = {0};
        jj =1; 
    
%         while jj <= numel(fileNames)
            for ii = 1:length(all_tets)
                % Collect spike data
                [spike_mat,spike_count] = extract_spikes(all_tets{ii},ii,timestamps,50); 
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
%     tetrode9 = cell2mat(new_tets(9,:)');
%     tetrode10 = cell2mat(new_tets(10,:)');
%     tetrode11 = cell2mat(new_tets(11,:)');
%     tetrode12 = cell2mat(new_tets(12,:)');    

    % Sum all spikes in trial
    num_spikes = sum(cell2mat(num_spikes),2);
    
    
    toc;
    %     save spikes and spike count from extract spikes output/ other variabes you want to save)
       
    
    


    