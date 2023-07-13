function collect_wfs(data)
    % loops through all the rows of spatData and uses
    % them to run dat2wfs on all the relevant .dats for all the trials. The
    % input -spatData
    % output - table with all the wf info collected for the
    % whole spatData - can include wf_means and max channel wfs - whatever
    % the waveformPCA code needs as input. - this is per rat and saves in a
    % folder called waveforms in the big computer 

    % Issue #1 - trialnames are missing but trial numbers can be used to 
    % run through the .dats in order 
    % Issue #2 - sleep spike times.  need to incorporate
    % concatenation of .dat plus find a way to get pre sleep filtering spike 
    % times - spike times from
    % unconcatenated flies can also be used if I put in a time offset -
    % just have to retreive the scan.mat for each animal 
    % Issue #3- metadata missing to choose mapping - will need to hardcode 

    load (data, 'spatData');

    % make index for single experiments on spatData 
    unique_dataset = unique(spatData.dataset);
    spatData_idx = [];
    for itS = 1: height(spatData)
        spatData_idx(itS) = find(strcmp(unique_dataset, spatData.dataset(itS)));
    end
    spatData_idx = spatData_idx'; 

    %make input variables for dat2wfs from spatData

    [read_dir, write_dir, mapping]= get_spatData_info(spatData);

    %check if any trials were concatenated and concatenate if
    %necessary.

    %if the trialname for sleep has more than 1
    %letter between the date and the underscore then exit matlab and
    %concatenate the sleep trials with names based on the letters in the gap - this needs to work for 2 or 3 trials.  

%     all_wfs = cell(size(spatData.trialNo)-1); %makes a cell array to be filled by wf info from all the files and all the trials -1 to remove sleep trial

    for itD = 1:(length(unique(spatData.dataset)))
        % need to loop through dataset variable to run dat2wfs per dataset
        % and the trial number one to run for all .dats in the folder 
        data_idx = find(spatData_idx == itD, 1,'first');
        data_idx_2 = find(spatData_idx == itD, 1,'last');
        exp_wfs = cell(length(data_idx:data_idx_2),size(spatData.trialNo,2)); %taking out -1 here so it picks up all 6 trials
        for trial_it = 1: (size(spatData.trialNo,2)) % -1 here temporarily removes the sleep trial there's one other joined file that might cause an issue which I changed manually  
            %insert sleep trial concatenation check loop here 
            if strcmp(string(spatData.env(data_idx,trial_it)), 'sleep')
                %extract sleep trial name
                sleepTrialName = spatData.trialName{data_idx,trial_it};
                sleepTrialCode = extractBefore(sleepTrialName, '_');
                sleepTrialDate = extractBefore(sleepTrialCode, 7); 
                sleepTrialLetters = extractAfter(sleepTrialCode, 6); %its always 6 because thats the data format
                if length(sleepTrialLetters) > 1 %if the sleep trial letters are longer than 1 concatenate the appropriate files and come back to running the code. 
                    % Run terminal command to concatenate multiple .dat files
                    dats_to_combine = [];
                    for letter_it = 1:length(sleepTrialLetters)
                        curr_sleep_dat = [read_dir{data_idx},sleepTrialDate, sleepTrialLetters(letter_it), '_sleepHP.dat'];
                        dats_to_combine = [dats_to_combine; curr_sleep_dat];
                    end
                    dats_to_combine = join(cellstr(dats_to_combine), ' ');
                    dats_to_combine = dats_to_combine{1}
                    combined_dat = [read_dir{data_idx} sleepTrialCode '_sleepHP.dat'];
                    command = ['cat ', dats_to_combine, ' > ' combined_dat '; exit'];
                    status = system(command);                 
                    if status == 0
                        disp('Concatenation completed successfully.');
                        % Continue with the rest of your MATLAB code
                    else
                        disp('Concatenation failed.');
                        % Handle the error or exit the program
                    end
                    %after concatenating the trial the trial iterator needs
                    %to become the right one for the concatenated trial 
                    trial_num =  trial_it + 1;% 
                     display(['two sleep trials, trial it: ', string(trial_it), 'trial_num: ' string(trial_num)])
                else 
                    trial_num =  trial_it; %this isn't working 
                    display(['one sleep trial, trial it: ', string(trial_it), 'trial_num: ' string(trial_num)])
                end 
            else
                trial_num = trial_it; 
                 display(['not a sleep trial, trial it: ', string(trial_it), 'trial_num: ' string(trial_num)])
            end

            [waveforms] = dat2wfs(read_dir{data_idx}, mapping{data_idx}, trial_num, trial_it, data, write_dir{data_idx},data_idx, data_idx_2); % run dat2wfs making wf. mats. not sure anymore if this should save a .mat anymore or just ouptut the wf_mat to be used here                    
            exp_wfs(:,trial_it) = waveforms; %save all wfs per experiment 
        end 
        save(fullfile('/data/isabella/probe_data/waveforms/', [spatData.dataset{data_idx} '_wfs.mat']), 'exp_wfs') %saves wfs per experiment 
    end
    
end

function [read_dir, write_dir, mapping]= get_spatData_info (spatData)
    % make all the input variables for dat2wfs from spatData 

    read_dir = cell(size(spatData.dataset,1),1);
    write_dir = cell(size(spatData.dataset,1),1);
    mapping = cell(size(spatData.dataset,1),1); 
    
    for itCl = 1:height(spatData) %iterate through cells in spatData

        %make read and write directories from spatData.dataset - hc
        %directories on big computer 
        
        dataset = spatData.dataset{itCl};
        animal = spatData.animal{itCl};
        date = extractBetween(dataset, '_', '_');
        read_dir{itCl}= ['/data/isabella/probe_data/' animal '/' date{1} '_1/'];
        write_dir{itCl} = ['/data/isabella/probe_data/' animal '/' date{1} '_1/writeDir/'];

        %load the correct mapping file - hc but should make a metadata
        %matrix per rat 
        if strcmp(animal, 'r1099')
            mapping{itCl} = 'map_single_same_plugbottom_final.mat';
        elseif strcmp(animal, 'r1143') || strcmp(animal, 'r1142')|| strcmp(animal, 'r1100')
            mapping{itCl} = 'map_single_opp_plugtop_final.mat';
        else 
            mapping{itCl} = 'map_multi_opp_plugtop_final.mat';
        end
    end

end