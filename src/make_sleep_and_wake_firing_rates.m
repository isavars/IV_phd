
function make_sleep_and_wake_firing_rates(data, clusters)

%load spatial data
load(data, 'spatData')

%load clusters 
load (clusters, 'PCA2_clusters','DG_ExCluster', 'CA3_ExCluster')

%make clusters 
granule_cluster = DG_ExCluster(PCA2_clusters == 2); %its 2 in this particular case (also for tetrode_adults_saved)
mossy_cluster = DG_ExCluster(PCA2_clusters == 1);
CA3_cluster = CA3_ExCluster; 

%make sleep variables 

%make indexes for sleep and wake trials 

maxWakeMeanRate= zeros(size(spatData,1),1);
sleep_idx = zeros(size(spatData,1),1);
awakeMeanRate = zeros(size(spatData,1),1);
for itCl = 1: height(spatData)
    sleep_trials = strcmp(string(spatData.env(itCl,:)),'sleep');
    sleep_idx_temp = find(sleep_trials);
    if size(sleep_idx_temp,2) > 1 
        sleep_idx(itCl) = sleep_idx_temp(2); %dealing with trial with more than one sleep
    else 
        sleep_idx(itCl) = sleep_idx_temp;
    end
    nov_trials = strcmp(string(spatData.env(itCl,:)),'nov');
    fam_trials = strcmp(string(spatData.env(itCl,:)),'fam');
    wake_trials = nov_trials + fam_trials; 
    %datasets have different numbers of wake trials 
    wake_idx_temp = find(wake_trials);
    maxWakeMeanRate(itCl) = nanmax(spatData.meanRate(itCl,wake_idx_temp));
    awakeMeanRate(itCl) = nanmean(spatData.meanRate(itCl,wake_idx_temp));
end 

%wake 
wakeMeanRate = maxWakeMeanRate; %maximum mean rate in wake trials 

%NREM - adapted for shifting sleep trial position 
NREMmeanRate = [];
for ii = 1:height(spatData)   
    curr_sleep_idx = sleep_idx(ii);
    NREMmeanRate = [NREMmeanRate ; spatData.meanRate(ii, curr_sleep_idx)]; %select sleep trial
end 


%REM
REMmeanRate = spatData.remMeanRate; 

%loop through groups and make plots for each
groups = {granule_cluster, mossy_cluster, CA3_cluster};
group_names = {'Granule Cells', 'Mossy Cells', 'CA3 Cells'};

%loop over spatData and assign excitatory cells to age bins

%gather age data from cellInfo 
cellInfo = getCellInfo(spatData);
age = cellInfo(:,3);

%separate the data into age bins 

ageBins = [17 20; 21 32]; %21 24; 

    for it_grp = 1:length(groups)
        original_group = groups{it_grp};
%         original_group = group;
        for itAge = 1: size(ageBins,1)            
            indAge = age>=ageBins(itAge,1) & age<=ageBins(itAge,2); % index for current age bin  
            %get indices for new clusters - shorten cell cluster to age bin
            %only. 
            ages_indexes_in_spatData = [];
            for it_indAge =1: length(indAge) 
                if true(indAge(it_indAge))
                    ages_indexes_in_spatData = [ages_indexes_in_spatData;it_indAge];
                else
                end
            end
            Newgroup = ismember(original_group,ages_indexes_in_spatData);
            group = original_group(Newgroup);
    
            % Calculate the mean value for "Wake" feature
            mean_values = [
                mean(wakeMeanRate(group)),    % Wake
                mean(NREMmeanRate(group)),    % NREM
                mean(REMmeanRate(group))      % REM
            ];
            
            % Calculate the standard deviation for each feature
            std_values = [
                std(wakeMeanRate(group)),    % Wake
                std(NREMmeanRate(group)),    % NREM
                std(REMmeanRate(group))      % REM
            ];
            
            % Define the names of the features
            feature_names = {'Wake', 'NREM', 'REM'};
            
            % Create the bar plot
            figure;
            bar(mean_values);
            hold on;
            
            % Add error bars
            errorbar(mean_values, std_values, '.');
            
            % Set the labels for the x-axis and y-axis
            xticks(1:numel(feature_names));
            xticklabels(feature_names);
            ylabel('firing rate (Hz)');

            %make Agebin titles 

            if itAge == 1
                AgeBin = 'Pre-wean';
            else 
                AgeBin = 'Post-wean';
            end
               
            
            % Set the title for the plot
            title(['Sleep Sate Firing Rates: ' group_names{it_grp} ' ' AgeBin]);
            
            hold off;

            %reset variables
%             group = original_group;
        end
    end
end     