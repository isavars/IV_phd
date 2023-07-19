function make_firing_property_plots(data,electrodes,cell_clusters) 
%this funciton makes plots showing different firig properties over
%age for the full population of cells (minus INS). will need DG_Ex_Cluster
%from class cells. boxplots for - mean firing rate during wake, sleep,
%sleep wake ratio, burst index. 

    %load spatData
    load (data, 'spatData');
    %load useful parts from spatData
    meanRate = spatData.meanRate;
    awakeMeanRate = nanmean(meanRate(:,1:5),2);
    burstIndex = spatData.burstIndex;
    dataset_spat = spatData.dataset;

    
    %get excitatory cell clusters 

    load(cell_clusters, 'PCA2_clusters', 'DG_ExCluster','CA3_ExCluster')     
    %make clusters from PCA2_clusters keeping old naming convention for
    %convenience in running the old code. 
    mossy_cluster =[];
    granule_cluster =[]; %keeping og naming convention for a second to see if this code runs 
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 2
            granule_cluster = [granule_cluster;DG_ExCluster(ii)]; %for WF1_AMR_BI_CCC_new - 2 is gc and 1 is mc 
        elseif PCA2_clusters(ii) == 1
            mossy_cluster = [mossy_cluster;DG_ExCluster(ii)];% for WF1_AMR_BI_DS2_new - 2 is gc and 1 is mc 
        end
    end

    pyramidal_cluster = CA3_ExCluster;

    cluster = {granule_cluster, mossy_cluster, pyramidal_cluster}; 
    clustername = {'granule', 'mossy', 'CA3'}; %to me used as tile change also depending on clusters plotted.         

%     %need to get this outside the loop - also this is the length of the
%     DG ex cluster - not spat data so the indexing needs adjusting 
% 
%     [~, ~, co_recorded_cells_capped] = class_cells(data,electrodes,cell_clusters);

    %loop over clusters and make sets of plots for each 
    for itC = 1:length(cluster)
        %make the features to be plotted 
        awakeMeanRate_clu = awakeMeanRate(cluster{itC});
        sleepMeanRate_clu = meanRate (cluster{itC},end);
        rateChange_clu = awakeMeanRate_clu ./ sleepMeanRate_clu;
        burstIndex_clu = burstIndex (cluster{itC});
%         co_recorded_cells_capped_clu = co_recorded_cells_capped(cluster{itC});
    
        %loop over spatData and assign excitatory cells to age bins
    
        %gather age data from cellInfo 
        cellInfo = getCellInfo(spatData);
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
        data_age_bin1 = [ awakeMeanRate_clu(bin_indices == 1), sleepMeanRate_clu(bin_indices == 1), rateChange_clu(bin_indices == 1), burstIndex_clu(bin_indices == 1)];
        data_age_bin2 = [ awakeMeanRate_clu(bin_indices == 2), sleepMeanRate_clu(bin_indices == 2), rateChange_clu(bin_indices == 2), burstIndex_clu(bin_indices == 2)];
%         data_age_bin3 = [ awakeMeanRate_clu(bin_indices == 3), sleepMeanRate_clu(bin_indices == 3), rateChange_clu(bin_indices == 3), burstIndex_clu(bin_indices == 3)];
        
        % Create labels for the title
        labels = {'awakeMeanRate', 'sleepMeanRate', 'rateChange', 'burstIndex'};
        
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
%     
%         corecorded_age_bin1 = co_recorded_cells_capped_clu(bin_indices == 1);
%         corecorded_age_bin2 = co_recorded_cells_capped_clu(bin_indices == 2);
%         corecorded_age_bin3 = co_recorded_cells_capped_clu(bin_indices == 3);
%     
%         figure;
%         bar([mean(corecorded_age_bin1), mean(corecorded_age_bin2), mean(corecorded_age_bin3)]);
%         hold on;
%         errorbar([mean(corecorded_age_bin1), mean(corecorded_age_bin2), mean(corecorded_age_bin3)], [std(corecorded_age_bin1), std(corecorded_age_bin2), std(corecorded_age_bin3)], '.');
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


