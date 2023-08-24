function plot_class_PCA_features_over_DS2( data, clusters, option)
    %this function needs to plot the features from the class PCA for a dataset
    %over ds2 orientation for that dataset - it can be boxplots but ideally
    %violin plots - the data i want to display is ds2 orientation on the x and
    %each feature on the y - the clusters will be side by side for each ds2
    %orientation. proportion of cells in each oreintaaiton should be counted
    %also. I want to plot post weaner only
    
    %load spatdata
    load(data,'spatData')   

    %load clusters 
    load(clusters, 'PCA2_clusters', 'DG_ExCluster') 
    

    %trim spatData to excitatory only 
    spatData = spatData(DG_ExCluster,:);

    %filter data by age - everything going in this funciton has to be excitatory
    %cells only 
    [age_filtered_spatData] = make_age_filtered_data(spatData, spatData);
    [age_filtered_PCA2_clusters] = make_age_filtered_data(spatData, PCA2_clusters);

    %make features 
    [burstIndex, awakeMeanRate, TP_Lat] = make_features(age_filtered_spatData);

    %obtain features from saved files 
    if option ==1 
        %load data from probes 
        load('probes_class_PCA_features.mat', 'DS2_orientations', 'co_recorded_shank_capped','wfPC1')
        [co_recorded_shank_capped] = make_age_filtered_data(spatData, co_recorded_shank_capped);
        [wfPC1] = make_age_filtered_data(spatData, wfPC1);
        features = {burstIndex, awakeMeanRate, wfPC1, co_recorded_shank_capped}; 
        feature_names = {'Burst Index','Mean Firing Rate (wake)','Wf-PCA PC1', 'Co-recorded Cells'};
    elseif option ==2 
        %load from tets 
        load ('tets_class_PCA_features.mat', 'DS2_orientations','ratio_silent_or_active_per_tet')
        [ratio_silent_or_active_per_tet] = make_age_filtered_data(spatData, ratio_silent_or_active_per_tet);
        features = {burstIndex, awakeMeanRate, TP_Lat, ratio_silent_or_active_per_tet};
        feature_names = {'Burst Index', 'Mean Firing Rate (wake)', 'Trough to Peak Latency', 'Ratio of Silent to Active'};
    end 
    
    %get DS2 orientations (new method)
    DS2_orientations = get_DS2_labels_by_distance_from_inversion();
    
    %get ds2 orientaitons
    [age_filtered_DS2_orientations] = make_age_filtered_data(spatData, DS2_orientations);

    %gather indexes for granule and mossy 
    granule_idx = find(age_filtered_PCA2_clusters == 1);
    mossy_idx = find(age_filtered_PCA2_clusters == 2);

   
    %loop through all features and make plots 
    for jj = 1:length (features)
        feature = features{jj};
        feature_name = feature_names{jj};
        %make features to plot   
        feature_GC_H = []; %granule cells in the Hillus
        feature_GC_G = []; %granule cells in the Granule cell layer
        feature_MC_H = []; %mossy cells in the Hillus
        feature_MC_G = []; %mossy cells in the Granule cell layer
        
        for ii = 1:length(feature)
            if find(ii == granule_idx)
                if age_filtered_DS2_orientations(ii) == 1 %pointing up
                    feature_GC_H = [feature_GC_H; feature(ii)]; 
                elseif  age_filtered_DS2_orientations(ii) == 2 || age_filtered_DS2_orientations(ii) == 3 %pointing down or inverting
                    feature_GC_G = [feature_GC_G; feature(ii)];
                end                 
            elseif find(ii == mossy_idx)
                if age_filtered_DS2_orientations(ii) == 1 %pointing up
                    feature_MC_H = [feature_MC_H; feature(ii)]; 
                elseif  age_filtered_DS2_orientations(ii) == 2 || age_filtered_DS2_orientations(ii) == 3 %pointing down or inverting
                    feature_MC_G = [feature_MC_G; feature(ii)];
                end 
            end 
        end 
    
        datasets = {feature_GC_G, feature_MC_G, feature_GC_H, feature_MC_H};
          
        %make boxplots and swarmcharts 
        make_plots(datasets,feature_name)
        %make half violins 
        %make_labia(datasets,feature_name)
    end


end 
function [burstIndex, awakeMeanRate, TP_latency] = make_features(spatData)
        %get data from maximum trial
        STs = spatData.SpkTs; 
        %burst index for clustering (Knierim method) 
        burstIndex = []; 
        for itC = 1: height(spatData)
            burstIndex = [burstIndex;(sum(diff(STs{itC}) <= 0.006))/(length(diff(STs{itC})))];% 0.008s for 
        end 
        %make awake mean rate  
        awakeMeanRate = zeros(size(spatData,1),1);
        wake_idx = cell(size(spatData,1),1);
        for itCl = 1: height(spatData)
            nov_trials = strcmp(string(spatData.env(itCl,:)),'nov');
            fam_trials = strcmp(string(spatData.env(itCl,:)),'fam');
            wake_trials = nov_trials + fam_trials; 
            %datasets have different numbers of wake trials 
            wake_idx_temp = find(wake_trials);
            wake_idx{itCl} = wake_idx_temp;
            awakeMeanRate(itCl) = nanmean(spatData.meanRate(itCl,wake_idx_temp));
        end 
        %make TP_lat
        for itSp = 1: length(spatData.nSpks) 
            [~, maxSpksPos] = nanmax(spatData.nSpks(itSp,:)); 
            TP_latency(itSp,:) = spatData.TP_latency(itSp, maxSpksPos);
        end
end 
function [age_filtered_data] = make_age_filtered_data(spatData, data)
    %get ages to make age bins
    cellInfo = getCellInfo(spatData);
    age = cellInfo(:,3);

    postwean_data_idx = [];
    for ii = 1: height(age)
        if age(ii) >= 21 
            postwean_data_idx = [postwean_data_idx; ii];
        end
    end 
    age_filtered_data = data(postwean_data_idx,:);

end
function make_labia(datasets,feature_name)
        
        % Set parameters
        boxwidth = 0.7; % Width of boxes
        linewidth = 2; % Line thickness for boxes and violins
        linecolor = 'k'; % Outline color
        palette = parula(numel(datasets)); % Color palette for violins
        
        % Function to calculate box-plot stats
        stat_fxn = @(x) [prctile(x,25), median(x), prctile(x,75), prctile(x,5), prctile(x,95)];
        
        % Prepare figure and subplots
        figure('color','w');
        plts = arrayfun(@(x) subplot(1, numel(datasets), x, 'nextplot', 'add'), 1:numel(datasets));
        
        % Generate violins
        [f, xi] = deal({}, {});
        for i = 1:numel(datasets)
            x = datasets{i};

            % Title for each subplot
            titles = {[feature_name ':GCL, Granule'], 'GCL, Mossy', 'Hilus, Granule', 'Hilus, Mossy'};
            title(plts(i), titles{i});
    
            % Skip if dataset is empty
            if isempty(x)
                continue;
            end
%             %trying out scaling 
%             f{i} = f{i} * length(x);

            [f{i}, xi{i}] = ksdensity(x, linspace(prctile(x, 1), prctile(x, 99), 100));
            patch(0 - [f{i}, zeros(1, numel(xi{i}), 1), 0], [xi{i}, fliplr(xi{i}), xi{i}(1)], palette(i, :), 'parent', plts(i));
        end
        
%         % Calculate the maximum extent for symmetry
%         maxval_list = zeros(1, numel(f));
%         for i = 1:numel(f)
%             if isnumeric(f{i}) && ~isempty(f{i})
%                 maxval_list(i) = max(f{i}(:));
%             else
%                 maxval_list(i) = -inf;  % or some other default value
%             end
%         end
%         maxval = max(maxval_list);
% 
% 
%         
%         % Generate box plots
%         for i = 1:numel(datasets)
%             x = datasets{i};
%             
%             % Skip if dataset is empty
%             if isempty(x)
%                 continue;
%             end
%             
%             stats = stat_fxn(x);
%             line([maxval/2, maxval/2], [stats(4), stats(5)], 'parent', plts(i));
%             patch(rescale([0, maxval, maxval, 0, 0], maxval * (1 - boxwidth), maxval * boxwidth), ...
%                 [stats(1), stats(1), stats(3), stats(3), stats(1)], palette(i, :), 'parent', plts(i));
%             line([maxval * (1 - boxwidth), maxval * boxwidth], [stats(2), stats(2)], 'parent', plts(i));    
%         end
%         
%         % Adjust properties for visualization
%         lines = findobj(gcf, 'Type', 'Line');
%         arrayfun(@(line) set(line, 'LineWidth', linewidth, 'Color', linecolor), lines);
%         
%         patches = findobj(gcf, 'Type', 'Patch');
%         arrayfun(@(patch) set(patch, 'LineWidth', linewidth, 'EdgeColor', linecolor), patches);
%         
%         [x0, x1] = deal(-maxval, maxval);
%         arrayfun(@(x) set(x, 'XLim', [x0, x1], 'xlabel', [], 'ylabel', [], 'xcolor', 'w', 'ycolor', 'w'), plts);

end 
function make_plots(datasets,feature_name)
    figure;
    % Define positions and labels
    positions = [2, 4, 7, 9];
    labels = {'GCL, Granule', 'GCL, Mossy', 'Hilus, Granule', 'Hilus, Mossy'}; 
    
    % Plot boxplots at specified positions
    hold on;
    for i = 1:4
        boxplot(datasets{i}, 'positions', positions(i), 'widths', 0.5);
    end
    
    % Adjust the x-axis properties
    set(gca, 'XTick', positions, 'XTickLabel', labels);
    title(feature_name);
    ylabel('Measurement Value');
    xlim([0, 11]);

    %make stacked histograms 

    figure;

    % Define labels
    labels = {'GCL, Granule', 'GCL, Mossy', 'Hilus, Granule', 'Hilus, Mossy'};
    
    % Define colors for each histogram for visibility
    colors = lines(numel(datasets));
    
    % Determine global x and y limits for consistent axis scaling
    globalXMax = -Inf;
    globalXMin = Inf;
    globalYMax = -Inf;
    
    for i = 1:numel(datasets)
        if ~isempty(datasets{i})
            globalXMax = max(globalXMax, max(datasets{i}));
            globalXMin = min(globalXMin, min(datasets{i}));
            
            counts = histcounts(datasets{i});
            globalYMax = max(globalYMax, max(counts));
        end
    end
    
    % Plot the first two datasets on the top subplot
    subplot(2, 1, 1); 
    hold on;
    for i = 1:2
        if ~isempty(datasets{i})
            histogram(datasets{i}, 'FaceColor', colors(i, :), 'FaceAlpha', 0.5);
        end
    end
    title('GCL Groups');
    ylabel('Frequency');
    legend(labels(1:2), 'Location', 'Best');
    xlim([globalXMin, globalXMax]);
    ylim([0, globalYMax]);
    
    % Plot the second two datasets on the bottom subplot
    subplot(2, 1, 2); 
    hold on;
    for i = 3:4
        if ~isempty(datasets{i})
            histogram(datasets{i}, 'FaceColor', colors(i, :), 'FaceAlpha', 0.5);
        end
    end
    title('Hilus Groups');
    ylabel('Frequency');
    xlabel('Measurement Value');
    legend(labels(3:4), 'Location', 'Best');
    xlim([globalXMin, globalXMax]);
    ylim([0, globalYMax]);
    
    sgtitle(feature_name); % Overall title for the entire figure


%     %make swarmchart
%     figure;
%     % Define positions and labels
%     positions = [1, 2, 4, 5];
%     labels = {'GCL, Granule', 'GCL, Mossy', 'Hilus, Granule', 'Hilus, Mossy'};
%     
%     % Combine datasets into one array and create a grouping variable
%     dataCombined = vertcat(datasets{:});
%     groups = arrayfun(@(x) x * ones(size(datasets{x})), 1:4, 'UniformOutput', false);
%     groupCombined = vertcat(groups{:});
%     
%     % Plot swarmchart
%     hold on;
%     swarmchart(groupCombined, dataCombined,5, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b','SizeData', 30); % 0.4 is the Jitter value, adjust as needed
%   
%     % Adjust the x-axis properties
%     set(gca, 'XTick', positions, 'XTickLabel', labels, 'YScale', 'log');
%     title(feature_name);
%     ylabel('Measurement Value');
%     xlim([0, 6]);
end