function make_firing_property_plots(data,electrodes) 
%this funciton needs to make plots showing different firig properties over
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

    %load elePos -might not be needed 
    load(electrodes, 'elePos')
    %load useful parts from elePos
    age = elePos.age;
    ds2_max_amplitude = nanmax(abs(elePos.DS2_mean_amplitude), [], 2);%make 1 max ds2_amplitude value per row of elePos

    %plots od DS2 amplitude over age 

    % Calculate mean and range per age
    [age_unique, ~, age_idx] = unique(age);
    mean_ds2_max_amplitude = splitapply(@nanmean, ds2_max_amplitude, age_idx);
    range_ds2_max_amplitude = splitapply(@(x) max(x) - min(x), ds2_max_amplitude, age_idx);
    
    % Plot mean with error bars
    figure;
    scatter(age_unique, mean_ds2_max_amplitude, 'filled', 'g');
    hold on;
    errorbar(age_unique, mean_ds2_max_amplitude, range_ds2_max_amplitude, 'LineStyle', 'none', 'Color', 'k');
    hold off;
    
    % Adjust x-axis range
    xMin = min(age_unique) - 1;
    xMax = max(age_unique) + 1;
    xlim([xMin, xMax]);
    
    % Add labels and title
    xlabel('Age');
    ylabel('Mean DS2 Max Amplitude');
    title('DS2 Max Amplitude over Age');

    
    %get excitatory cells from putative DG only 
    [~,~, DG_ExCluster, ~] = class_cells(data,electrodes);
    
    %make the features to be plotted 
    awakeMeanRate = awakeMeanRate(DG_ExCluster);
    sleepMeanRate = meanRate (DG_ExCluster,end);
    rateChange = awakeMeanRate ./ sleepMeanRate;
    burstIndex = burstIndex (DG_ExCluster);

    %loop over spatData and asign excitatory cells to age bins

    %gather age data from cellInfo 
    cellInfo = getCellInfo(spatData);
    age = cellInfo(:,3);

    %separate the data into age bins 
    
    age_bins = [17 20; 21 24; 25 31];
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
    
    
    % Define the data for each age bin
    data_age_bin1 = [awakeMeanRate(bin_indices == 1), sleepMeanRate(bin_indices == 1), rateChange(bin_indices == 1), burstIndex(bin_indices == 1)];
    data_age_bin2 = [awakeMeanRate(bin_indices == 2), sleepMeanRate(bin_indices == 2), rateChange(bin_indices == 2), burstIndex(bin_indices == 2)];
    data_age_bin3 = [awakeMeanRate(bin_indices == 3), sleepMeanRate(bin_indices == 3), rateChange(bin_indices == 3), burstIndex(bin_indices == 3)];
    
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
            data = {data_age_bin1(:, i), data_age_bin2(:, i), data_age_bin3(:, i)}; % };%
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
        title(sprintf('Boxplot for %s', labels{i}));
%         ylabel('Value');
        xticks(1:num_groups);
        xticklabels(age_bin_labels(1:num_groups));
%         xtickangle(45);
    end





% Define the variables we want to test
properties = {awakeMeanRate, sleepMeanRate, rateChange, burstIndex};
property_names = {'awakeMeanRate', 'sleepMeanRate', 'rateChange', 'burstIndex'};

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
    [p, ~, stats] = kruskalwallis(data_vector, group_vector, 'off');
    fprintf('Kruskal-Wallis test for %s: p = %.4f\n', property_names{i}, p);
    
    figure;
    % Post hoc test (e.g., pairwise comparisons)
    c = multcompare(stats, 'CType', 'dunn-sidak');  % Pairwise comparison with Dunn-Sidak correction
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
