
function make_sleep_and_wake_firing_rates(data, clusters)

%load spatial data
load(data, 'spatData')

%load clusters 
load (clusters, 'PCA2_clusters','DG_ExCluster', 'CA3_ExCluster')

%make clusters 
granule_cluster = DG_ExCluster(PCA2_clusters == 2); %its 2 in this particular case 
mossy_cluster = DG_ExCluster(PCA2_clusters == 1);
CA3_cluster = CA3_ExCluster; 

%make sleep variables 
%wake 
wakeMeanRate = spatData.meanRate(:,1:5); %maximum mean rate in wake trials 
for it_mr = 1: height (spatData)
    wakeMeanRate(it_mr,:) = max(wakeMeanRate(it_mr,:),[],'omitnan');
end
wakeMeanRate = wakeMeanRate(:,1); 

%NREM
NREMmeanRate = spatData.meanRate(:,6); %sleep trial 

%REM
REMmeanRate = spatData.remMeanRate; 

%loop through groups and make plots for each
groups = {granule_cluster, mossy_cluster, CA3_cluster};
group_names = {'Granule Cells', 'Mossy Cells', 'CA3 Cells'};

    for it_grp = 1:length(groups)

        % Calculate the mean value for "Wake" feature
        mean_values = [
            mean(wakeMeanRate(groups{it_grp})),    % Wake
            mean(NREMmeanRate(groups{it_grp})),    % NREM
            mean(REMmeanRate(groups{it_grp}))      % REM
        ];
        
        % Calculate the standard deviation for each feature
        std_values = [
            std(wakeMeanRate(groups{it_grp})),    % Wake
            std(NREMmeanRate(groups{it_grp})),    % NREM
            std(REMmeanRate(groups{it_grp}))      % REM
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
        
        % Set the title for the plot
        title(['Sleep Sate Firing Rates: ' group_names{it_grp}]);
        
        hold off;

    end

end     