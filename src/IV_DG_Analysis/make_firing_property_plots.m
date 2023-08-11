function make_firing_property_plots(data,electrodes,cell_clusters) 
%this funciton makes plots showing different firig properties over
%age for the full population of cells (minus INS). will need DG_Ex_Cluster
%from class cells. boxplots for - mean firing rate during wake, sleep,
%sleep wake ratio, burst index. 

    %load spatData
    load (data, 'spatData');
    %load useful parts from spatData
    meanRate = spatData.meanRate;
%     awakeMeanRate = nanmean(meanRate(:,1:5),2);
    burstIndex = spatData.burstIndex;


    %make indexes for sleep and wake trials 

    maxWakeMeanRate= zeros(size(spatData,1),1);
    sleep_idx = zeros(size(spatData,1),1);
    awakeMeanRate = zeros(size(spatData,1),1);
    wake_idx = cell(size(spatData,1),1);
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
        wake_idx{itCl} = wake_idx_temp;
        maxWakeMeanRate(itCl) = nanmax(spatData.meanRate(itCl,wake_idx_temp));
        awakeMeanRate(itCl) = nanmean(spatData.meanRate(itCl,wake_idx_temp));

    end 
     


    %gather age data from cellInfo 
    cellInfo = getCellInfo(spatData);


    
    %get excitatory cell clusters 

    load(cell_clusters, 'PCA2_clusters', 'DG_ExCluster','CA3_ExCluster')     
    %make clusters from PCA2_clusters keeping old naming convention for
    %convenience in running the old code. 
    mossy_cluster =[];
    granule_cluster =[]; %keeping og naming convention for a second to see if this code runs 
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 1
            granule_cluster = [granule_cluster;DG_ExCluster(ii)]; %for WF1_AMR_BI_CCC_new - 2 is gc and 1 is mc 
        elseif PCA2_clusters(ii) == 2
            mossy_cluster = [mossy_cluster;DG_ExCluster(ii)];% for WF1_AMR_BI_DS2_new - 2 is gc and 1 is mc 
        end
    end

    pyramidal_cluster = CA3_ExCluster;

%     cluster = {granule_cluster, mossy_cluster, pyramidal_cluster}; 
%     clustername = {'granule', 'mossy', 'CA3'}; %to me used as title change also depending on clusters plotted.  

    cluster = {granule_cluster, mossy_cluster, DG_ExCluster}; 
    clustername = {'granule', 'mossy', 'Pooled Data'}; %to me used as title change also depending on clusters plotted.  

%     %need to get this outside the loop - also this is the length of the
%     DG ex cluster - not spat data so the indexing needs adjusting 
% 
%     [~, ~, co_recorded_cells_capped] = class_cells(data,electrodes,cell_clusters);

    [SVAR_G,SVAR_M,SVAR_C ]= make_silent_vs_active_per_tet(spatData, cellInfo, cluster, sleep_idx, wake_idx );


    %loop over clusters and make sets of plots for each 
    for itC = 1:length(cluster)
        %make the features to be plotted 
        awakeMeanRate_clu = awakeMeanRate(cluster{itC});
        sleepMeanRate_clu = meanRate (cluster{itC},sleep_idx(cluster{itC}));
        rateChange_clu = awakeMeanRate_clu ./ sleepMeanRate_clu;
        burstIndex_clu = burstIndex (cluster{itC});
%         co_recorded_cells_capped_clu = co_recorded_cells_capped(cluster{itC});
        if itC ==1
            SVAR_clu = SVAR_G; %silent to active ratio granule 
        elseif itC ==2
            SVAR_clu = SVAR_M; %silent to active ratio mosssy
        elseif itC ==3
            SVAR_clu = SVAR_C; %silent to active ratio ca3 cells 
        end 

        %loop over spatData and assign excitatory cells to age bins
        
        %get age from cell info
        age = cellInfo(:,3);
        age = age(cluster{itC});
   
        %separate the data into age bins 
        
        age_bins = [17 20; 21 31]; %21 24; 
        bin_indices = zeros(size(age, 1), 1);
        
        for itA = 1:size(age, 1)
            current_age = age(itA);  
            for itBins = 1:size(age_bins, 1)
                if current_age >= age_bins(itBins, 1) && current_age <= age_bins(itBins, 2)
                    bin_indices(itA) = itBins;
                    break;
                end
            end
        end
        
        %make box plots for continous variables 
        
        % Define the data for each age bin
        data_age_bin1 = [ awakeMeanRate_clu(bin_indices == 1), SVAR_clu(bin_indices == 1), rateChange_clu(bin_indices == 1), burstIndex_clu(bin_indices == 1)];
        data_age_bin2 = [ awakeMeanRate_clu(bin_indices == 2), SVAR_clu(bin_indices == 2), rateChange_clu(bin_indices == 2), burstIndex_clu(bin_indices == 2)];
%         data_age_bin3 = [ awakeMeanRate_clu(bin_indices == 3), sleepMeanRate_clu(bin_indices == 3), rateChange_clu(bin_indices == 3), burstIndex_clu(bin_indices == 3)];
        
        % Create labels for the title
        labels = {'awakeMeanRate', 'slient_to_acive_ratio', 'rateChange', 'burstIndex'};
        
        % Create labels for age bins
        age_bin_labels = cell(1, size(age_bins, 1));
        for i = 1:size(age_bins, 1)
            age_bin_labels{i} = sprintf('P %d to P %d', age_bins(i, 1), age_bins(i, 2));
        end
        
        % Plot separate boxplots for each variable and age bin
        figure;
        for i = 1:4
            %subplot(2, 2, i);
            figure;
            if i <= 3
                data = {data_age_bin1(:, i), data_age_bin2(:, i)};%, data_age_bin3(:, i)}; % 
                num_groups = numel(data);
                positions = 1:num_groups;
                for j = 1:num_groups
                    boxplot(data{j}, 'Positions', positions(j), 'Whisker', Inf);
                    hold on;
                end
                set(gca, 'YScale', 'log'); % Set y-axis to log scale
            else
                data = {data_age_bin1(:, i), data_age_bin2(:, i)};%, data_age_bin3(:, i)};
                num_groups = numel(data);
                positions = 1:num_groups;
                for j = 1:num_groups
                    boxplot(data{j}, 'Positions', positions(j), 'Whisker', Inf);
                    hold on;
                end
            end
            hold off;
            title(['Boxplot for', labels{i}, ': ', clustername{itC}]);
    %         ylabel('Value');
            xticks(1:num_groups);
            xticklabels(age_bin_labels(1:num_groups));
    %         xtickangle(45);
        end
%     
%         %make bar plots for discrete variables 
    
%         corecorded_age_bin1 = co_recorded_cells_capped_clu(bin_indices == 1);
%         corecorded_age_bin2 = co_recorded_cells_capped_clu(bin_indices == 2);
% %         corecorded_age_bin3 = co_recorded_cells_capped_clu(bin_indices == 3);
%     
%         figure;
%         bar([mean(corecorded_age_bin1), mean(corecorded_age_bin2)]);%, mean(corecorded_age_bin3)]);
%         hold on;
%         errorbar([mean(corecorded_age_bin1), mean(corecorded_age_bin2)], [std(corecorded_age_bin1), std(corecorded_age_bin2)], '.');
%         %errorbar([mean(corecorded_age_bin1), mean(corecorded_age_bin2), mean(corecorded_age_bin3)], [std(corecorded_age_bin1), std(corecorded_age_bin2), std(corecorded_age_bin3)], '.');
%         set(gca, 'XTickLabel', {'Pre-wean', 'Post-wean'});
%         ylabel('Mean co-recorded cells');
%         title('Mean co-recorded cells per age');
        
    %         normality tests 
    %         % Create histograms for each firing property
    %         for i = 1:4
    %             figure;
    %             histogram([data_age_bin1(:, i); data_age_bin2(:, i); data_age_bin3(:, i)], 'Normalization', 'probability');
    %             title(sprintf('Histogram for %s', labels{i}));
    %             xlabel('Value');
    %             ylabel('Probability');
    %         end
        
        %run kruskal-wallis for non-parametric data swap out variables 
    
        % Define the variables we want to test
        properties = { awakeMeanRate_clu, sleepMeanRate_clu, rateChange_clu, burstIndex_clu};%, co_recorded_cells_capped_clu};
        property_names = {'awakeMeanRate', 'sleepMeanRate', 'rateChange', 'burstIndex'};%,'co_recorded cells'};
    
        KW_for_age_groups(age_bins, age, properties, property_names)
    end

end

function KW_for_age_groups(age_bins, age, properties, property_names)

    % Loop over each property
    for i = 1:numel(properties)
    
        % Initialize the data and groups for the Kruskal-Wallis test
        data_vector = [];
        group_vector = [];
    
        % Loop over each age bin
        for j = 1:size(age_bins, 1)
    
            % Find the indices of the observations in the current age bin
            bin_indices = (age >= age_bins(j, 1)) & (age <= age_bins(j, 2));
            
            % Append the observations to the data vector
            data_vector = [data_vector; properties{i}(bin_indices)];
            
            % Append the group identifiers to the group vector
            group_vector = [group_vector; repmat(j, sum(bin_indices), 1)];
        end
    
        % Perform the Kruskal-Wallis test
        [p, ~, stats] = kruskalwallis(data_vector, group_vector, 'off')
        fprintf('Kruskal-Wallis test for %s: p = %.4f\n', property_names{i}, p);
        
        figure;
        % Post hoc test (e.g., pairwise comparisons)
%         c = multcompare(stats, 'CType', 'dunn-sidak');  % Pairwise comparison with Dunn-Sidak correction
    end  

end 

function [SVAR_G,SVAR_M,SVAR_C ]= make_silent_vs_active_per_tet(spatData, cellInfo, cluster, sleep_idx, wake_idx )
%this needs to read in cells from amb cluster and provide say if they are active in sleep only or at least one wake env 
% then take the average mean firing rate of all cells on that tet and rank
% them from low to high - introdicting cuttoffs where the ratio of active
% in wake vs only active in sleep flips - that can be considered the hillus
% and then when it filps back it can be considered CA3 - go per rat and do
% all depths then see if it aligns with DS2 - can only spend time on this
% if you find a moment. 

    %itterate over clusters and find if the cells are active in sleep only or
    %in run trials for each cluster 
    for itClu = 1:length(cluster)
        curr_cluster = cluster{itClu};
        spatData_temp = spatData(curr_cluster,:);    
        silent_or_active = zeros(height(spatData_temp),1);
        for itC = 1: height(spatData_temp)
            if all(spatData.nSpks(itC,wake_idx{itC}) < 75)
                silent_or_active(itC) = 0; %active in sleep only
            else 
                silent_or_active(itC) = 1; %active in wake trial
            end
        end
        %get the mean rate per tetrode 
        meanRate_per_tet = zeros(height(spatData_temp),1);
        ratio_silent_or_active_per_tet = zeros(height(spatData_temp),1);
        cellInfo_temp = cellInfo(curr_cluster,:);
        for jj=1:length(cellInfo_temp)
            ID = cellInfo_temp(jj, 1);
            tet = cellInfo_temp(jj, 2);
            age = cellInfo_temp(jj, 3);
            dates = cellInfo_temp(jj, 5);      
            % Find indices of rows with the same ID, tet, age, and date
            matchingIndices = find(cellInfo_temp(:,1) == ID & cellInfo_temp(:,2) == tet & cellInfo_temp(:,3) == age & cellInfo_temp(:,5) == dates);
            %get mean rate per tet 
            meanRate_per_tet(jj) = nanmean(spatData.meanRate(matchingIndices,sleep_idx(jj)));
            ratio_silent_or_active_per_tet(jj) = sum(silent_or_active(matchingIndices) ==0)/sum(silent_or_active(matchingIndices) ==1);
            if ratio_silent_or_active_per_tet(jj) == inf
                ratio_silent_or_active_per_tet(jj) = 1;
            end
        end  
        if itClu == 1 
            SVAR_G = ratio_silent_or_active_per_tet;
        elseif itClu == 2
            SVAR_M = ratio_silent_or_active_per_tet;
        elseif itClu == 3
            SVAR_C = ratio_silent_or_active_per_tet;
        end
    end
end


    %normality tests 
%         % Create histograms for each firing property
%         for i = 1:4
%             figure;
%             histogram([data_age_bin1(:, i); data_age_bin2(:, i); data_age_bin3(:, i)], 'Normalization', 'probability');
%             title(sprintf('Histogram for %s', labels{i}));
%             xlabel('Value');
%             ylabel('Probability');
%         end
%         
%         
%         % Create P-P plots for each firing property
%         for i = 1:4
%             figure;
%             data_combined = [data_age_bin1(:, i); data_age_bin2(:, i); data_age_bin3(:, i)];
%             normplot(data_combined);
%             title(sprintf('P-P Plot for %s', labels{i}));
%             xlabel('Expected Cumulative Probability');
%             ylabel('Observed Cumulative Probability');
%         end


