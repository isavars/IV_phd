function [new_elePos]= makeElePos(rat_information, readDir)
    %make elepos from DS2_info files. one row is one DS2 info file so per rat
    %per age. Age and trialname should be in the table. 
    % TO DO 1) make sure this works with a tetrode file DS2_info

    %load all the DS2_info files in readDir
    fileName = dir(strcat(readDir, '/*.mat'));  %finds all the .mat in the folder
   
    %make a new table: elePos 
    varNames = {'rat_ID', 'probe_type', 'hist_labels', 'dataset', 'age', 'DS2_channels','DS2_spike_times' 'DS2_labels', 'DS2_max_amplitude', 'DS2_mean_amplitude','DS2_rate','DS2_peak_to_trough_amplitude','DS2_slope'};  
    varTypes = {'string', 'double', 'string', 'string', 'double', 'double', 'cell','string','double', 'double', 'double','double', 'double'};
    new_elePos = table('Size', [height(fileName), length(varTypes)], 'VariableNames', varNames, 'VariableTypes', varTypes);
    new_elePos.hist_labels = repmat("", height(fileName), 4);
    new_elePos.DS2_channels = repmat(nan, height(fileName), 32);
    new_elePos.DS2_labels = repmat("",height(fileName), 32);
    new_elePos.DS2_max_amplitude = repmat(nan, height(fileName), 32);
    new_elePos.DS2_mean_amplitude = repmat(nan, height(fileName), 32);
    new_elePos.DS2_peak_to_trough_amplitude = repmat(nan, height(fileName), 32);
    new_elePos.DS2_slope = repmat(nan, height(fileName), 32);  
    new_elePos.DS2_spike_times = cell(height(fileName),1); 
    %loop through all the rat IDs in elePos and add data for first few rows
    % Get the list of files in the directory
    fileList = dir(fullfile(readDir, '*.mat'));
    load(rat_information, 'rat_info');
    % Iterate over the files
    for i = 1:numel(fileList)
        % Get the current filename
        filePath = fullfile(readDir, fileName(i).name);
        load(filePath, 'DS2_info')
        % Extract the relevant information from the loaded file and assing
        % to new elepos
        parts = split(fileName(i).name, '_');
        rat_ID = parts{1};
        dataset_name = extractBefore(fileName(i).name, '_DS2_info.mat');
        age = str2double(extractBetween( fileName(i).name, 'P', '_DS2_info.mat'));
        new_elePos.dataset(i) = dataset_name;
        new_elePos.age(i) = age;
        % Find rows with matching rat_id
        % Convert rat_id column in ElePos to strings
        elePos_rat_id = rat_info.rat_id;
        elePos_rat_id = cellstr(elePos_rat_id);  % convert string array to cell array of character vectors
        elePos_rat_id = cellfun(@(x) strrep(x,' ',''), elePos_rat_id, 'UniformOutput', false);
        elePos_rat_id = string(elePos_rat_id);

        % Convert rat_id_to_search to string if it's not already
        rat_ID = string(rat_ID);
        matching_rows = elePos_rat_id == rat_ID; 
        new_elePos.rat_ID(i) = rat_ID;
        % Retrieve the associated probe_type and hist_labels and Assign the extracted values to the new_elePos table
        matching_probe_type = rat_info.probe_type(matching_rows);
        matching_hist_labels = rat_info.hist_labels(matching_rows, :);
        new_elePos.probe_type(i) = matching_probe_type;
        new_elePos.hist_labels(i,:) = matching_hist_labels;
        new_elePos.DS2_channels(i,:) = DS2_info.channel';
        new_elePos.DS2_labels(i,:) = DS2_info.DS2_labels';
        new_elePos.DS2_max_amplitude(i,:) = DS2_info.DS2_max_amplitude';
        new_elePos.DS2_mean_amplitude(i,:) = DS2_info.DS2_mean_amplitude';
        new_elePos.DS2_peak_to_trough_amplitude(i,:) = DS2_info.DS2_peak_to_trough_amplitude';
        new_elePos.DS2_slope(i,:) = DS2_info.DS2_slope';
        new_elePos.DS2_rate(i,:) = DS2_info.DS2_rate(1);
        new_elePos.DS2_spike_times{i} = DS2_info.DS2_spike_times(1);
%         for jj = 1:32
%             new_elePos.DS2_spike_times{jj} = DS2_info.DS2_spike_times;
%         end
    end
end 