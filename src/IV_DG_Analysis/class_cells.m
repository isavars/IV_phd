function [DS2_amplitudes,PCA2_clusters, DG_ExCluster, InCluster] = class_cells(data,electrodes)
    %this function needs to read in cells from spatData and electrode
    %position table (which has DS2 info and other postion estimate methods)
    %and provide cell classification ids (the length of spatData) for 
    %interneuron, granule cell, mossy cell and CA3 pyramidal cell. 

    % TO DO - Add linear probe info and distance from inversion info to
    %         DS2_info.
    %       - fix DS2_amplitude to be a mean of values around the peak
    %       instead of the peak 

    %load spatial Data 
    load (data, 'spatData');
    %load useful parts from spatData
    meanRate = spatData.meanRate;
    awakeMeanRate = nanmean(meanRate(:,1:5),2);
    burstIndex = spatData.burstIndex;
    animal = spatData.animal;
    dataset = spatData.dataset;
    wf_means = spatData.wf_means;
    waveforms = spatData.waveforms; %waveforms is created from the mean wf with the maximum amplitude - need to keep channel id for it to use DS2 labels
    max_wf_chan = spatData.max_wf_channel; %this is per trial needs to be one value 
    nSpks = spatData.nSpks;
    SpkTs = spatData.SpkTs;
    TP_latency = mean(spatData.TP_latency,2);
    %load electrode Positions 
    load (electrodes, 'elePos');
    rat_id = elePos.rat_ID;% cellstr(elePos.rat_id); the output from new_elePos is slightly different to elePos
    dataset_elePos = (elePos.dataset);
    hist_labels = elePos.hist_labels;

    %chose trial with most spikes to use wfs, spike times, mean firing rate, max wf channel from it 
    for itSp = 1: length (nSpks) 
        [~, maxSpksPos] = max(nSpks(itSp,1:5)); % all wake trials (hc change somewhere)
        WFs (itSp,:) = wf_means(itSp, maxSpksPos); %gets best wf from wf means
        max_wf_chan(itSp,:) = max_wf_chan(itSp,maxSpksPos); 
        STs (itSp,:) = SpkTs(itSp, maxSpksPos); %gets the spiketimes from the same trial to calculate the burst index
%         awakeMeanRate (itSp,:) = meanRate(itSp, maxSpksPos); %gets the wake mean firing rate for the most active wake trial (doesn this make sense?)
    end
    
    %get tetrode to shank and channel mapping to be used by cluster feature
    %functions 
    [tetShankChan] = makeTetShankChan(spatData,max_wf_chan);


    %step 1 - filter by gross histology postion CA3 vs DG 
        % reads elePos values for hist_labels per tetrode and if a cell is
        % in the DG it gets labeled DG_cluster and if its in CA3
        % CA3_cluster. 
    
     DG_cluster = [];
     CA3_cluster = []; 

     for it_cells = 1: height(spatData)
         for it_ep = 1: height(elePos)
             if strcmp(animal(it_cells),rat_id(it_ep)) 
                 hist_label = hist_labels(it_ep,tetShankChan(it_cells,2));
                 if strcmp(hist_label, "DG")
                    DG_cluster = [DG_cluster;it_cells];
                 elseif strcmp(hist_label, "CA3")
                    CA3_cluster = [CA3_cluster;it_cells]; 
                 end
             end
         end
     end  
    % this is a branching point for DG and CA3 - the CA3 then needs to be
    % sorted into interneurons and CA3 pyramidal cells also - should It all
    % be in one big loop that iterates over rows of spatData? 


    %step 2 - filter by excitatory vs inhibitory using trough-to-peak
    %measure, mean rate, and maybe burst index from Knierim. 

        %make the inhibitory and excitatory clusters for CA3 and DG subgroups
        %-mantaining row numbers from spatData
        InCluster = [];
        ExCluster = [];
        for itWF = 1: length (WFs)
            if TP_latency(itWF) < 0.425 && awakeMeanRate(itWF) > 1.2 
                InCluster = [InCluster;itWF]; 
            else 
                ExCluster = [ExCluster;itWF]; 
            end            
        end 
        % temporarily making a DG_Ex_cluster
        DG_ExCluster = DG_cluster(ismember(DG_cluster, ExCluster));

    %step 3 - call in position relative to DS2 (specifically DS2 amplitude
    %for the ones that have it and more informed histology lables for the
    %ones that dont) 
        %loop through spatData(DG_ExCluster) and find the DS2 amplitudes 
        % of the chosen cells per channel with the largest wf.  

        DS2_amplitudes = []; %trying this out with mean amplitude 
        for it_DG_Ex = DG_ExCluster'
            for it_ep = 1: height(elePos)
                if strcmp(dataset(it_DG_Ex),dataset_elePos(it_ep))
                    DS2_amplitude = elePos.DS2_amplitude(it_ep,tetShankChan(it_DG_Ex,3)); %third column of tetShankChan is for channel
                    DS2_amplitudes = [DS2_amplitudes;DS2_amplitude];
                end
            end 
        end

    %step 4 - create wf_PCA features (call waveformPCA, run as a
    %subfunciton) need to know what outputs from the pca are the right ones
    %to use here. 

        [wfPC1,wfPC2] = waveformPCA(DG_ExCluster,waveforms);

    %step 5 - create sleep vs wake firing rate (rateChange) - needs to be
    %indexed within the cluster for DG or CA3 and withought Interneurons -
    %need to make a DG_exCluster and a CA3_exCluster - or do this
    %differently with one big loop 

        awakeMeanRate = awakeMeanRate(DG_ExCluster);
        sleepMeanRate = meanRate (DG_ExCluster,end);
        rateChange = awakeMeanRate ./ sleepMeanRate;
    
    %step 6 - make "slope" from Knierim group (slope of best fit line 
    % through normalized, sorted waveform peaks of the four tetrode wires)

        slope = zeros(length(DG_ExCluster), 1);        
        for i = 1:length(DG_ExCluster)            
            waveform = WFs{DG_ExCluster(i)};% Select the waveform for the current index            
            peaks = max(waveform, [], 1);% Determine the peak values for each tetrode channel         
            sorted_peaks = sort(peaks, 'descend');% Sort the peak values in ascending order  
            normalized_peaks = sorted_peaks./max(sorted_peaks); %(sorted_peaks - mean(sorted_peaks)) / std(sorted_peaks);% Normalize the sorted peak values      
            x = 1:length(normalized_peaks);% Calculate the slope of the differences using linear regression
            coeffs = polyfit(x, normalized_peaks, 1);
            slope(i) = abs(coeffs(1));
        end

    %step 7 - make co-recorded cells -for a given cell how many other cells
    %were recorded on that shank or tetrode - overlap tetrodes might
    %complicate this. use tetShankChan. 
    %make a measure thats co-recorded so that any cell recorded on the same
    %shank as 4+ cells gets grouped into one 

        co_recorded_cells = zeros(size(tetShankChan, 1), 2);
        for i = 1:size(tetShankChan, 1)
            tetrode = tetShankChan(i, 1);
            shank = tetShankChan(i, 2);
            dataset_id = dataset{i};% Get the dataset identifier for the current cell            
            % Find the rows in tetShankChan with the same tetrode number and dataset identifier
            same_tetrode_rows = (tetShankChan(:, 1) == tetrode) & strcmp(dataset_id, dataset);            
            % Count the number of co-recorded cells on the same tetrode (excluding the current cell)
            co_recorded_tetrode_count = sum(same_tetrode_rows) - 1;            
            % Find the rows in tetShankChan with the same shank number and dataset identifier
            same_shank_rows = (tetShankChan(:, 2) == shank) & strcmp(dataset_id, dataset);            
            % Count the number of co-recorded cells on the same shank (excluding the current cell)
            co_recorded_shank_count = sum(same_shank_rows) - 1;            
            % Store the co-recorded counts for the current cell
            co_recorded_cells(i, 1) = co_recorded_tetrode_count;
            co_recorded_cells(i, 2) = co_recorded_shank_count;
        end
        
        co_recorded_shank = co_recorded_cells(DG_ExCluster,2);

        co_recorded_shank_capped = zeros(size(co_recorded_shank,2));
        for ii = 1:length(co_recorded_shank)
            if co_recorded_shank(ii) >= 4 
                co_recorded_shank_capped(ii) = 4;
            else 
                co_recorded_shank_capped(ii) = co_recorded_shank(ii);
            end
        end 
        co_recorded_shank_capped = co_recorded_shank_capped';

    %step 8 - run a second PCA on the chosen features -output is PC1 to be
    %used in clustering and histogram showing the distribution. 
        
        %burst index for clustering (Knierim method) 
        burstIndex = []; 
        for it_DE = DG_ExCluster'
            burstIndex = [burstIndex;(sum(diff(STs{it_DE}) <= 0.008))/(length(diff(STs{it_DE})))];% 0.008s produced the best bimodality 
        end 
%         burstIndex = burstIndex(DG_ExCluster);

        data = [wfPC1, rateChange, burstIndex,co_recorded_shank_capped];        % Combine the variables into a matrix
        [PC1, PC2]= class_PCA(data);        % run PCA
        
    % step 9 -run k means with PCs from second PCA and other features 
        %trying a k-means on different features 
        cluster_data = [PC1,PC2];%         cluster_data = [wfPC1,PC1];
        %normalize variances 
        cluster_data = cluster_data./nanstd(cluster_data,0,1);                 
        [PCA2_clusters]= kmeans_clustering(cluster_data); %try different features in here 

%         % Perform k-means clustering on the four features in data (that
%         made the best second PCA PC1. 
%         k = 2;  % Number of clusters
%         idx = kmeans(data, k);
%         PCA2_clusters = idx;
  

    %testing plots     

    % Define the cluster labels
    cluster1 = PCA2_clusters == 1;
    cluster2 = PCA2_clusters == 2;
    
    % Plot the histogram of the first principal component
    figure;
    hold on;
    h1 = histogram(PC1(cluster1), 'NumBins', 12, 'FaceColor', 'blue');
    h2 = histogram(PC1(cluster2), 'NumBins', 12, 'FaceColor', 'red');
    hold off;
    
    xlabel('First Principal Component');
    ylabel('Frequency');
    title('Histogram of the First Principal Component');
    
    % Create a legend for the clusters
    legend([h1 h2], 'Cluster 1', 'Cluster 2');

    
%     %k-means on awake mean rate and wfPC1 - this might still be coming up
%     with the best clusters ... 
%         cluster_data = [wfPC1,awakeMeanRate];
%         %normalize variances 
%         cluster_data = cluster_data./nanstd(cluster_data,0,1);
%         % Run k-means clustering
%         k = 2;  % Set the number of clusters to 2
%         [idx, centroids] = kmeans(cluster_data, k);    
%         % creating the clustering variable 
%         wfPC_rate_clusters =idx;
%         % Plotting the clusters
%         figure;
%         gscatter(data(:, 1), data(:, 2), idx);
%         hold on;
%         scatter(centroids(:, 1), centroids(:, 2), 100, 'k', 'filled');
%         legend('Cluster 1', 'Cluster 2', 'Centroids'); 
%         xlabel("wfPC1")
%         ylabel("wake rate")
%         set(gca,'Yscale','log')
    
%     figure;
%     gscatter(wfPC1, DS2_amplitudes, PCA2_clusters, 'bg', '.', 12);
%     xlabel("wfPC1")
%     ylabel("DS2 Amplitudes")
% 
%     figure;
%     gscatter(awakeMeanRate, DS2_amplitudes, PCA2_clusters, 'bg', '.', 12);
%     xlabel("firing rate")
%     set(gca,'Xscale','log')
%     ylabel("DS2 Amplitudes")
% 
%     figure;
%     gscatter(burstIndex, DS2_amplitudes, PCA2_clusters, 'bg', '.', 12);
%     xlabel("burstIndex")
%     ylabel("DS2 Amplitudes")
% 
%     figure;
%     hold all;
%     gscatter(co_recorded_shank, DS2_amplitudes, PCA2_clusters, 'bg', '.', 12);
%     xlabel("# co-recorded cells on the same shank")
%     ylabel("DS2 Amplitudes")
%        
%     figure;
%     hold all;
%     gscatter(wfPC1, awakeMeanRate, wfPC_clusters, 'bg', '.', 12);
%     xlabel("wfPC1")
%     ylabel("wake rate")
%     set(gca,'Yscale','log')
% 
%     figure;
%     hold all;
%     gscatter(burstIndex, awakeMeanRate, wfPC_rate_clusters, 'rg', '.', 12);
%     xlabel("burst index")
%     ylabel("wake rate")
%     set(gca,'Yscale','log')
% 
%     figure;
%     hold all;
%     gscatter(burstIndex, wfPC1, wfPC_rate_clusters, 'rg', '.', 12);
%     xlabel("burst index")
%     ylabel("wfPC1")
% 
%     figure;
%     hold all;
%     gscatter(co_recorded_shank, wfPC1, wfPC_rate_clusters, 'rg', '.', 12);
%     xlabel("# co-recorded cells on the same shank")
%     ylabel("wfPC1")
% 
%     figure;
%     hold all;
%     gscatter(burstIndex, awakeMeanRate, PCA2_clusters, 'bg', '.', 12);
%     xlabel("burst index")
%     ylabel("wake rate")
%     set(gca,'Yscale','log')
% 
%     figure;
%     hold all;
%     gscatter(burstIndex, wfPC1, PCA2_clusters, 'bg', '.', 12);
%     xlabel("burst index")
%     ylabel("wfPC1")
% 
%     figure;
%     hold all;
%     gscatter(co_recorded_shank, wfPC1, PCA2_clusters, 'bg', '.', 12);
%     xlabel("# co-recorded cells on the same shank")
%     ylabel("wfPC1")


%     figure;
%     hold all;
%     gscatter(slope, awakeMeanRate, wfPC_rate_clusters, 'bg', '.', 12);
%     xlabel("slope")
%     ylabel("wake rate")
%     set(gca,'Yscale','log')
% 
%     figure;
%     hold all;
%     gscatter(slope, wfPC1, wfPC_rate_clusters, 'bg', '.', 12);
%     xlabel("slope")
%     ylabel("wfPC1")


end 

function [tetShankChan] = makeTetShankChan(spatData,max_wf_chan)
    %make tetrode labels from CellIDs in spatData - can probably write
        %something else for this 
        cellInfo = getCellInfo(spatData);
        tet = cellInfo(:,2);
        % create channel inex per tetrode 'tet_index' which is a tetrode by
        % channel identity array
        tet_index = zeros(8, 4); % create an 8x4 matrix of zeros
        tet_index(1:2:end, :) = repmat(1:8:25, 4, 1).' + repmat(0:2:6, 4, 1); % fill odd-numbered rows with odd numbers
        tet_index(2:2:end, :) = repmat(2:8:26, 4, 1).' + repmat(0:2:6, 4, 1); % fill even-numbered rows with even numbers
        overlap = repmat (0:8:24,4,1).' + repmat (5:8,4,1);%overlap tetrodes are the bottom 4 contacts which are closer to each other 
        tet_index = [tet_index; overlap]; % + 1 overlap per octrode
        
        tetShankChan = zeros(length(tet),3);
        % this needs to be per probe type
        for it_tet = 1: length(tet)
            tetShankChan(it_tet,1) = tet(it_tet);
            tetShankChan(it_tet,3) = tet_index(tet(it_tet),max_wf_chan(it_tet));
            if tet(it_tet) == 1 || tet(it_tet) == 2 || tet(it_tet) == 9
                tetShankChan(it_tet,2) = 1;
            elseif tet(it_tet) == 3 || tet(it_tet) == 4 || tet(it_tet) == 10
                tetShankChan(it_tet,2) = 2;
            elseif tet(it_tet) == 5 || tet(it_tet) == 6 || tet(it_tet) == 11
                tetShankChan(it_tet,2) = 3;
            elseif tet(it_tet) == 7 || tet(it_tet) == 8 || tet(it_tet) == 12
                tetShankChan(it_tet,2) = 4;
            end
        end
end 

function [wfPC1,wfPC2] = waveformPCA(DG_ExCluster,waveforms)

    ex_waveforms = waveforms(DG_ExCluster);
    pca_data = [];
    for itEx = 1: length(ex_waveforms)
         if ~isnan(cell2mat(ex_waveforms(itEx))) % && length(cell2mat(ex_waveforms(itEx))) == 97 %changed from 50 don't know why this is here 
             pca_data = [pca_data; interp1(1:97, cell2mat(ex_waveforms(itEx)),1:0.48:97,'spline')];%[pca_data; interp1(1:50, cell2mat(ex_waveforms(itEx)),1:0.48:50,'spline')]
         end
    end
    % Normalise so each WF peak is 1.
    pca_data      = pca_data ./ max(pca_data,[],2);
    
    % Shift WFs so that peaks aligned. - output starts at aligned peak
    % instead of data begining 
    [~,pkInd]  = max(pca_data,[],2);
    firstPkInd = min(pkInd);
    nSampsWF   = size(pca_data,2);
    padWF      = zeros(size(pca_data));
    maxShift = max(pkInd - firstPkInd); 
    for itWF=1:size(pca_data,1)
        WFShift = pkInd(itWF) - firstPkInd;
        padWF(itWF, (1+maxShift-WFShift):(nSampsWF-WFShift+maxShift)) = pca_data(itWF, :);
    end
    
    pca_data = padWF( : , firstPkInd:nSampsWF );%check ig this line works 
    
    % create a wave width at half height feature?

    % second derrivative 
    diff_pca_data = [];
    for itPCA = 1: size(pca_data,1)
        diff_pca_data = [diff_pca_data; diff(pca_data(itPCA,maxShift:end),2)]; %inserted 1:80 for the .8 ms time window % better time window - maxShift:end
    end
    [coeff,score] =  pca(diff_pca_data(:,1:160));
    w_PC1 = coeff(:,1);
    w_PC2 = coeff(:,2);

    
    %try just plotting sections of the second derrivative to see if you get
    %a bimodality 

    % weight of the PCA
    wfPC1 = score(:,1);
    wfPC2 = score(:,2);

    %k-means on the 1st two PCs of the wfPCA - not using for clustering 
    
%     data = [wfPC1,wfPC2];
%     %normalize variances 
%     data = data./nanstd(data,0,1);
%     % Run k-means clustering
%     k = 2;  % Set the number of clusters to 2
%     [idx, centroids] = kmeans(data, k);
%     
%     % creating the clustering variable 
%     wfPC_clusters =idx;
% 
%     % Plotting the clusters
%     figure;
%     gscatter(data(:, 1), data(:, 2), idx);
%     hold on;
%     scatter(centroids(:, 1), centroids(:, 2), 100, 'k', 'filled');
%     legend('Cluster 1', 'Cluster 2', 'Centroids');  

end
function [first_pc, second_pc] = class_PCA(data)
    % Remove rows with NaN values in any of the variables
    data(isinf(data)) = NaN;
    valid_rows = ~any(isnan(data), 2);
    data = data(valid_rows, :);

    % Apply log normalization to awakeMeanRate
    awakeMeanRate = data(:, 2);
    log_normalized_awakeMeanRate = log1p(normalize(awakeMeanRate));

    % Combine the log-normalized feature with the rest of the data
    transformed_data = [data(:, 1), log_normalized_awakeMeanRate, data(:, 3), data(:, 4)];

    % Normalize the data
    normalized_data = zscore(transformed_data);

    % Perform PCA on normalized data
    [coeff, score, ~, ~, explained] = pca(normalized_data);

    % Extract the first principal component
    first_pc = NaN(size(data, 1), 1);
    first_pc(valid_rows) = score(:, 1);
    
    % Extract the second principal component
    second_pc = NaN(size(data, 1), 1);
    second_pc(valid_rows) = score(:, 2);

    % Display the explained variance ratios
    explained_ratio = explained / sum(explained);
    disp('Explained Variance Ratios:');
    disp(explained_ratio);
end

function [PCA2_clusters]= kmeans_clustering(cluster_data)

        % Run k-means clustering
        k = 2;  % Set the number of clusters to 2
        [idx, centroids] = kmeans(cluster_data, k);    
        % creating the clustering variable 
        PCA2_clusters =idx;
        % Plotting the clusters
        figure;
        gscatter(cluster_data(:, 1), cluster_data(:, 2), idx);
        hold on;
        scatter(centroids(:, 1), centroids(:, 2), 100, 'k', 'filled');
        legend('Cluster 1', 'Cluster 2', 'Centroids'); 
        xlabel('feature 1')
        ylabel('feature 2')
        title('K-Means Clustering');
end



%%% Features I'm not using for now but might be useful to keep 

%         %make burst index from Senzai and Buzsaki 
%         burstIndex_Senzai = [];
%         for itSP = 1: length(SpkTs)
%             spike_AC = spk_crosscorr(cell2mat(STs(itSP)),'AC',0.001,0.3,900, 'norm', 'none'); 
%             spike = mean(spike_AC(304:306)); % maybe change this to 3-8ms
%             baseline= mean(spike_AC(501:601));
%             burstIndex_Senzai = [burstIndex_Senzai; spike/baseline];
%         end     

% testing plots 

%     figure;
%     hold all;
%     gscatter(wfPC1, DS2_amplitudes, wfPC_rate_clusters, 'br', '.', 12);
%     xlabel("wfPC1")
%     ylabel("DS2 Amplitudes")

% 
%     figure;
%     scatter (wfPC2,DS2_amplitudes)
%     xlabel("wfPC2")
%     ylabel("DS2 Amplitudes")
% 
%     figure;
% %     scatter (rateChange,DS2_amplitudes)
%     gscatter(rateChange, DS2_amplitudes, wfPC_clusters, 'br', '.', 12);
%     xlabel("firing rate change")
%     set(gca,'Xscale','log')
%     ylabel("DS2 Amplitudes")
% % 
%     figure;
% %     scatter (burstIndex,DS2_amplitudes)
%     gscatter(burstIndex, DS2_amplitudes, wfPC_clusters, 'br', '.', 12);
%     xlabel("burstIndex")
%     ylabel("DS2 Amplitudes")
%     
%     figure;
% %     scatter3(rateChange,DS2_amplitudes, wfPC1)
%     gscatter3(rateChange, DS2_amplitudes,wfPC1, wfPC_clusters, 'br', '.', 12);
%     xlabel("firing rate change")
%     ylabel("DS2 Amplitudes")
%     zlabel("wfPC1")

% %     Extract data for each cluster based on wfPC_clusters
%     cluster1_data = wfPC1(wfPC_clusters == 1);
%     cluster2_data = wfPC1(wfPC_clusters == 2);
% %     Plot overlayed histograms with transparent bars and full-color edges
%     figure;
%     hold all;
%     histogram(cluster1_data, 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'LineWidth', 2);
%     histogram(cluster2_data, 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2);    
%     xlabel('w-PC1');
%     ylabel('Cell Count');
%     title('Histogram of PC1');
%     legend('Cluster 1', 'Cluster 2');
% 
%     figure;
%     DS2_amplitudes_c1 = DS2_amplitudes(wfPC_clusters == 1);
%     DS2_amplitudes_c2 = DS2_amplitudes(wfPC_clusters == 2);
%     hold all;
%     histogram(DS2_amplitudes_c1, 'FaceAlpha', 0.5, 'EdgeColor', 'b', 'LineWidth', 2, 'NumBins',15);
%     histogram(DS2_amplitudes_c2, 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2,'NumBins',12); 
%     xlabel('DS2 Amplitudes');
%     ylabel('Frequency');
%     title('Histogram of DS2 Amplitudes');

%     figure()
%     histogram(wfPC1)
%     xlabel('w-PC1')
%     ylabel('Cell Count')
%     title('Histogram of PC1')
%     
%     figure()
%     histogram(wfPC2)
%     xlabel('w-PC2')
%     ylabel('Cell Count')
%     title('Histogram of PC2')

%     figure;
%     scatter(wfPC1,wfPC2)
%     xlabel("wfPC1")
%     ylabel("wfPC2")
% 
%     figure;
%     scatter (wfPC1,DS2_amplitudes)
%     xlabel("wfPC1")
%     ylabel("DS2 Amplitudes")

% function [first_pc, second_pc] = class_PCA(data)
%     % Remove rows with NaN values in any of the variables
%     data = data(~any(isnan(data), 2), :);
% 
%     % Apply log normalization to awakeMeanRate
%     awakeMeanRate = data(:, 2);
%     log_normalized_awakeMeanRate = log1p(normalize(awakeMeanRate));
% 
%     % Combine the log-normalized feature with the rest of the data
%     transformed_data = [data(:, 1), log_normalized_awakeMeanRate, data(:, 3), data(:, 4)];
% 
%     % Normalize the data
%     normalized_data = zscore(transformed_data);
% 
%     % Perform PCA on normalized data
%     [coeff, score, ~, ~, explained] = pca(normalized_data);
% 
%     % Extract the first principal component
%     first_pc = score(:, 1);
%     second_pc = score(:, 2);
% 
%     % Plot the histogram of the first principal component
%     figure;
%     histogram(first_pc,'NumBIns',20);
%     xlabel('First Principal Component');
%     ylabel('Frequency');
%     title('Histogram of the First Principal Component');
% 
%     % Display the explained variance ratios
%     explained_ratio = explained / sum(explained);
%     disp('Explained Variance Ratios:');
%     disp(explained_ratio);
% end