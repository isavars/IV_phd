function PCA2_clusters = class_cells_tetrodes(data,electrodes,cluster_filename)
    %this function needs to read in cells from spatData and electrode
    %position table (which has DS2 info and other postion estimate methods)
    %and provide cell classification ids (the length of spatData) for 
    %interneurons, granule cells, mossy cells and CA3 pyramidal cells. 
    %essentailly the same function as class_cells but with modifications
    %for tetrode data. 

    % TO DO:
    %       -make all the steps that arent subfuntions into subfunctions
    %       -check if wf-PCA with shorter wf length is working - need
    %       tetrode spatData for this. 
    %       -elepos needs to have channel number feature to assoicate ds2's
    %       with what tetrode they came from and be able to use the DS2
    %       orientation feature - wont be right as is. 


    %load spatial Data 
    load (data, 'spatData');

    %make indexes for sleep and wake trials 
    sleepMeanRate_all= zeros(size(spatData,1),1);
    sleep_idx = zeros(size(spatData,1),1);
    awakeMeanRate_all = zeros(size(spatData,1),1);
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
        awakeMeanRate_all(itCl) = nanmean(spatData.meanRate(itCl,wake_idx_temp));
        sleepMeanRate_all(itCl) = nanmean(spatData.peakRate(itCl,sleep_idx_temp));
    end 

    %remove outliers from the dataset - i belive i did a bad job cutting
    %cells so this is here - cells need to have a minimum firing rate during sleep 
    new_spatData = [];
    for ii = 1:height(spatData)
        if sleepMeanRate_all(ii) > 0.01 %&& spatData.burstIndex(ii,sleep_idx(ii)) > 0.00001%&& sleepMeanRate_all(ii) > awakeMeanRate_all(ii)%
            new_spatData = [new_spatData; spatData(ii,:)];
        end
    end 
    spatData = new_spatData;
    height(spatData)

    %load useful parts from spatData
    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    dataset = spatData.dataset;
    wf_means = spatData.wf_means;
    waveforms = spatData.waveforms; %waveforms is created from the mean wf with the maximum amplitude - need to keep channel id for it to use DS2 labels
    max_wf_chan = spatData.max_wf_channel; %this is per trial needs to be one value 
    nSpks = spatData.nSpks;
    SpkTs = spatData.SpkTs;
    %TP_latency = nanmean(spatData.TP_latency,2);
    %load electrode Positions -rename this to metadata or add the option
    %for elePos or a different metadata file to be used. 
    load (electrodes, 'elePos');
    rat_id = elePos.rat_ID;% cellstr(elePos.rat_id); the output from new_elePos is slightly different to elePos
    dataset_elePos = (elePos.dataset);
    hist_labels = elePos.hist_labels; 
    probe_type = elePos.probe_type;

    %chose trial with most spikes to use wfs, spike times, mean firing rate, max wf channel from it 
    for itSp = 1: length (nSpks) 
        [~, maxSpksPos] = nanmax(nSpks(itSp,:)); 
%         WFs (itSp,:) = wf_means(itSp, maxSpksPos); %gets best wf from wf means
        max_wf_chan(itSp,:) = max_wf_chan(itSp,maxSpksPos); 
        max_waveforms(itSp,:) = waveforms(itSp,maxSpksPos); 
        STs (itSp,:) = SpkTs(itSp, maxSpksPos); %gets the spiketimes from the same trial to calculate the burst index
        max_burstIndex(itSp,:) = burstIndex(itSp, maxSpksPos); 
        TP_latency(itSp,:) = spatData.TP_latency(itSp, maxSpksPos);
%         max_awakeMeanRate (itSp,:) = meanRate(itSp, maxSpksPos); %gets the wake mean firing rate for the most active wake trial (doesn this make sense?)
    end
       
    %make an equivalent of shank_channels is just tetrode the cell was
    %recorded on - is there an issue with 64 channels vs 32 anywhere? 
    cellInfo = getCellInfo(spatData);
    tets = cellInfo(:,2); %actually also need tetrodes from the channel in DS2_channels - so that DS2 signals are mapped to channels 

    %make indexes for sleep and wake trials 

    maxWakeMeanRate= zeros(size(spatData,1),1);
    sleep_idx = zeros(size(spatData,1),1);
    awakeMeanRate_all = zeros(size(spatData,1),1);
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
        awakeMeanRate_all(itCl) = nanmean(spatData.meanRate(itCl,wake_idx_temp));
        %make WFs from wake trials only
        [~, maxPos] = nanmax(nSpks(itCl,wake_idx_temp));
        WFs(itCl,:) = spatData.wf_means(itCl,maxPos); %wfs come form wake trials 
    end 
    
    %this should be remove rate differenced caused by devleopmental changes
    %in ps
    %awakeMeanRate_all = maxWakeMeanRate;


    %step 1 - filter by gross histology postion CA3 , DG or ambiguous
        % reads elePos values for hist_labels - hist labels has 4
        % labels but this doesnt matter for the tetrodes because there is
        % no mapping 

     DG_cluster = [];
     CA3_cluster = []; 
     AMB_cluster = []; 

     for it_cells = 1: height(spatData)
         for it_ep = 1: height(elePos)
             if strcmp(dataset(it_cells),dataset_elePos(it_ep)) 
                 hist_label = hist_labels(it_ep,1); %changed to 1 because all the labels are the same for a dataset 
                 if strcmp(hist_label, "DG")
                    DG_cluster = [DG_cluster;it_cells];
                 elseif strcmp(hist_label, "CA3")
                    CA3_cluster = [CA3_cluster;it_cells]; 
                 elseif strcmp(hist_label, "AMB") %adding ambiguous cells 
                    AMB_cluster = [AMB_cluster;it_cells];                 
                 end
             end
         end
     end  


    %step 2 - filter by excitatory vs inhibitory using trough-to-peak
    %measure, mean rate, and maybe burst index from Knierim.     
    [DG_ExCluster, CA3_ExCluster, AMB_ExCluster, InCluster1, InCluster2, low_narrow] = interneuron_filter(WFs,TP_latency,awakeMeanRate_all,max_burstIndex,DG_cluster,CA3_cluster, AMB_cluster);

    
    %filter ambiguous clusters for CA3 vs DG cells using rem sleep - so in
        %S&B they show GCs have a much higher firing rate in NREM than wake or
        %REM - but its significantly different bettween all three mossy cells
        %are almost the same for all three but a bit more active in sleep and
        %CA3 cells are more active in wake and nrem than rem 


        %need to plot this again with tetrode data and see if it makes
        %sense - another option is just using DS2 amplidude and frequency /
        %are he cells firing around DS2. 
    
        %use AMB_ExCluster - loop through and sort cells into DG vs CA3 based
        %on mean rate differences between wake, nrem and rem. 
    
        %make sleep variables 
        
        %wake 
        wakeMeanRate = maxWakeMeanRate; %spatData.meanRate(:,wake_idx); %maximum mean rate in wake trials      
        amb_wakeMeanRate = wakeMeanRate(AMB_ExCluster);%wake 
        
        %make nrem mean rate for only amb cluster with varying index for
        %sleep trial 
        amb_nremMeanRate = [];
        for ii = 1:height(AMB_ExCluster)   
            curr_sleep_idx = sleep_idx(ii);
            amb_nremMeanRate = [amb_nremMeanRate ; spatData.meanRate(ii, curr_sleep_idx)]; %select sleep trial
        end 
        %rem rate 
        amb_remMeanRate = spatData.remMeanRate(AMB_ExCluster); %rem 
     
%         for it_amb = 1: length(AMB_ExCluster)
%             if amb_remMeanRate(it_amb)  < amb_wakeMeanRate(it_amb) && amb_remMeanRate(it_amb) < amb_nremMeanRate(it_amb) %behaves like a CA3 cell in wake vs sleep
%                 %add index from AMB_ExCluster into CA3_ExCluster
%                 CA3_ExCluster = [CA3_ExCluster;AMB_ExCluster(it_amb)]; 
%             else 
%                 %add it to DG_ExCluster
%                 DG_ExCluster = [DG_ExCluster;AMB_ExCluster(it_amb)];
%             end
%         end  
%         %you will need to sort the clusters in ascending order after 
%         DG_ExCluster = sort(DG_ExCluster);
%         CA3_ExCluster = sort(CA3_ExCluster); 
        
        %atempt to sort amb clusters (new labels)
         [meanRate_per_tet,ratio_silent_or_active_per_tet]= make_silent_vs_active_per_tet(spatData, AMB_cluster,cellInfo, DG_ExCluster,amb_wakeMeanRate, sleep_idx, wake_idx );


    
    %step 3 - Position relative to DS2 (specifically DS2 direction
    %per shank as a discrete variable - for a given channel number produces
    %all the channel numbers on the same shank and decide if its pointing 
    % up or down through 3 possibe options for calculating it - check
    %literature?  

        option = 3; %find the right one and the right thresholds 
    
        DS2_orientations = get_DS2_orientations(spatData,DG_ExCluster,elePos, tets, option);
        DS2_orientations = DS2_orientations(DS2_orientations~=0);

        %filter amb cluster using DS2 orientation 



    %step 4 - create wf_PCA features (call waveformPCA, run as a
    %subfunciton) need to know what outputs from the pca are the right ones
    %to use here. 

        [wfPC1,wfPC2, pca_data, diff_pca_data] = waveformPCA(DG_ExCluster,max_waveforms); % i think this needs to be the max waveform 

    
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

    %step 6- create sleep vs wake firing rate (rateChange) - needs to be
    %indexed within the cluster for DG or CA3 and withought Interneurons -
    %need to make a DG_exCluster and a CA3_exCluster - or do this
    %differently with one big loop 

        awakeMeanRate = awakeMeanRate_all(DG_ExCluster);
        sleepMeanRate = sleepMeanRate_all(DG_ExCluster);
        rateChange = awakeMeanRate ./ sleepMeanRate;

        %rateChange can have inf values if dividing by zero in sleep and 0
        %if rate is 0 in wake 
        for ii = 1: length(rateChange)
            if rateChange(ii) == inf || rateChange(ii) == 0 %these should basically not be in the dataset but if they are they shouldn't be classified as mossy
                rateChange(ii) = 1; % if its 0 would be off during wake and on in sleep so should be GC 
            end 
        end

    %step 7 - make co-recorded cells -for a given cell how many other cells
    %were recorded on that shank or tetrode - overlap tetrodes might
    %complicate this. use tetShankChan. 
    %make a measure thats co-recorded so that any cell recorded on the same
    %shank as 4+ cells gets grouped into one 
        
        co_recorded_tet = cellInfo(:,4);
        co_recorded_tet_DG_all = co_recorded_tet(DG_ExCluster); %makes co_recorded_cells all cells 

        %remove interneurons from co_recorded cells to make co recorded
        %excitatory only like in Knierim (try this one on class_cells also)     
        %cut down spatData to Excitartory cells only and create index for
        %rows as you go

        excitatory_spatData = spatData(DG_ExCluster,:);     
        co_recorded_tet_DG_ex_cell_info = getCellInfo(excitatory_spatData); %makes co_recorded_cells excitatory only 
        co_recorded_tet_DG_ex = co_recorded_tet_DG_ex_cell_info(:,4);

        co_recorded_tet_DG_ex_capped = zeros(size(co_recorded_tet_DG_ex,1),1); 
            for ii = 1:length(co_recorded_tet_DG_ex)
                if co_recorded_tet_DG_ex(ii) >= 4 
                    co_recorded_tet_DG_ex_capped(ii) = 4;
                else 
                    co_recorded_tet_DG_ex_capped(ii) = co_recorded_tet_DG_ex(ii);
                end
            end 

        co_recorded_tet_capped = zeros(size(co_recorded_tet_DG_all,1),1); 
        for ii = 1:length(co_recorded_tet)
            if co_recorded_tet(ii) >= 4 
                co_recorded_tet_capped(ii) = 4;
            else 
                co_recorded_tet_capped(ii) = co_recorded_tet(ii);
            end
        end 
%         co_recorded_tet_capped = co_recorded_tet_capped'; %capped at 4 to reduce variance in the group provided by outliers 
        co_recorded_tet_capped = co_recorded_tet_capped(DG_ExCluster);
    %step 8 - run a second PCA on the chosen features -output is PC1 to be
    %used in clustering and histogram showing the distribution. 
        
        %burst index for clustering (Knierim method) 
        burstIndex = []; 
        for it_DE = DG_ExCluster'
            burstIndex = [burstIndex;(sum(diff(STs{it_DE}) <= 0.006))/(length(diff(STs{it_DE})))];% 0.008s produced the best bimodality 
        end 
%         burstIndex = burstIndex(DG_ExCluster);
        
        TP_lat = TP_latency(DG_ExCluster);
        %normalize data for PCA
        [normalized_TP_lat] = normalize_data_by_age(TP_lat, DG_ExCluster, spatData);
        [normalized_awakeMeanRate] = normalize_data_by_age(sleepMeanRate, DG_ExCluster, spatData);
        [normalized_sleepMeanRate] = normalize_data_by_age(awakeMeanRate, DG_ExCluster, spatData);
        [normalized_rateChange] = normalize_data_by_age(rateChange, DG_ExCluster, spatData);
        [normalized_wfPC1] = normalize_data_by_age(wfPC1, DG_ExCluster, spatData);
        [normalized_wfPC2] = normalize_data_by_age(wfPC2, DG_ExCluster, spatData);

        data = [ normalized_TP_lat , awakeMeanRate, burstIndex, ratio_silent_or_active_per_tet];% Final combo
        %data = [ normalized_wfPC1 , normalized_rateChange, normalized_wfPC2, DS2_orientations];%Buzsaki combo
        %data = [ slope , normalized_sleepMeanRate, burstIndex,co_recorded_tet_DG_ex];%Knierim combo
        [PC1, PC2]= class_PCA(data);        % run PCA
        
    % step 9 -run k means with PCs from second PCA and other features 
        %trying a k-means on different features 
        cluster_data = [PC1,PC2];%         cluster_data = [wfPC1,PC1];
        %normalize variances 
        cluster_data = cluster_data./nanstd(cluster_data,0,1);                 
        [PCA2_clusters]= kmeans_clustering(cluster_data); %try different features in here 
        %save clusters 
        save(cluster_filename,'PCA2_clusters','DG_ExCluster','CA3_ExCluster', 'InCluster1', 'InCluster2', 'low_narrow')
        

    %testing plots - this needs to be a sub function at the very bottom 

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

    %wf-pca plots 
    
    %plot of the normalized waveforms 
    figure;
    hold all;
    colors = {'b','g'};
    for ii = 1:length(DG_ExCluster)
        cluster = PCA2_clusters(ii); 
        plot(pca_data(ii,:),'Color', colors{cluster});
    end

    %plot of the second derrivative of the normalized waveforms 
    figure;
    hold all;
    colors = {'b','g'};
    for ii = 1:length(DG_ExCluster)
        cluster = PCA2_clusters(ii); 
        plot(diff_pca_data(ii,:),'Color', colors{cluster});
    end
    hold off;

    %interneuron filter summary figure
    % Initialize an array with the same size as your data.
    groups = categorical(repmat("", size(TP_latency)));    
    % Set the groups based on the index arrays.
    groups(DG_ExCluster) = 'DG';
    groups(AMB_ExCluster) = 'Unclassified';
    groups(InCluster1) = 'InCluster1';
    groups(InCluster2) = 'InCluster2';    
    % Manually specify the order of the groups.
    groups = reordercats(groups, {'InCluster1', 'InCluster2', 'DG', 'Unclassified'});    
    % Define colors.
    colors = [
        0 0 1;   % dark blue for 'InCluster1'
        0.4 0.4 1; % light blue for 'InCluster2'
        0 1 0;   % green for 'DG'
        1 0 0;   % red for 'CA3'
    ];    

    % Create the scatter plot.
    figure;
    gscatter(TP_latency, awakeMeanRate_all, groups, colors);
    xlabel("Trough to Peak Latency (ms)",'FontSize', 16)
    ylabel("Mean Firing Rate (Hz)",'FontSize', 16)
    set(gca,'Yscale','log')
    % Add a legend.
    legend('InCluster1', 'InCluster2', 'DG', 'Unclassified');    

    % Create the scatter plot.
    figure;
    gscatter(max_burstIndex, awakeMeanRate_all, groups, colors);
    xlabel("Burst index",'FontSize', 16)
    ylabel("Mean Firing Rate (Hz)",'FontSize', 16)
    set(gca,'Xscale','log')
    set(gca,'Yscale','log')
    % Add a legend.
    legend('InCluster1', 'InCluster2', 'DG', 'Unclassified');     
    
    
    %plots for clustering features overlayed on DS2 orientations 
    figure;
    gscatter(wfPC1, jitter(DS2_orientations), PCA2_clusters, 'bg', '.', 12);
    xlabel("wfPC1",'FontSize', 16)
    ylabel("DS2 orientations",'FontSize', 16)
    set(gca, 'YTick', [1,2,3,4], 'YTickLabel', {'up', 'inverting', 'down', 'no DS2'},'FontSize', 16);

    figure;
    gscatter(awakeMeanRate, jitter(DS2_orientations), PCA2_clusters, 'bg', '.', 12);
    xlabel("firing rate")
    set(gca,'Xscale','log')
    ylabel("DS2 orientations")
    set(gca, 'YTick', [1,2,3,4], 'YTickLabel', {'up', 'inverting', 'down', 'no DS2'});

    figure;
    gscatter(burstIndex, jitter(DS2_orientations), PCA2_clusters, 'bg', '.', 12);
    xlabel("burstIndex")
    set(gca,'Xscale','log')
    ylabel("DS2 orientations")
    set(gca, 'YTick', [1,2,3,4], 'YTickLabel', {'up', 'inverting', 'down', 'no DS2'});

    figure;
    gscatter(slope, jitter(DS2_orientations), PCA2_clusters, 'bg', '.', 12);
    xlabel("slope")
    ylabel("DS2 orientations")
    set(gca, 'YTick', [1,2,3,4], 'YTickLabel', {'up', 'inverting', 'down', 'no DS2'});

    figure;
    hold all;
    gscatter(jitter(co_recorded_tet_DG_ex), jitter(DS2_orientations), PCA2_clusters, 'bg', '.', 12);
    xlabel("# co-recorded cells on the same tetrode")
    ylabel("DS2 orientations")
    set(gca, 'YTick', [1,2,3,4], 'YTickLabel', {'up', 'inverting', 'down', 'no DS2'});


end  

function [DG_ExCluster, CA3_ExCluster, AMB_ExCluster, InCluster1, InCluster2, low_narrow] = interneuron_filter(WFs,TP_latency, awakeMeanRate,max_burstIndex, DG_cluster,CA3_cluster, AMB_cluster)
        %make the inhibitory and excitatory clusters for CA3 and DG subgroups
        %-mantaining row numbers from spatData
        InCluster1 = []; %narrow wf INs 
        InCluster2 = []; %wide wf INs
        low_narrow = []; %confunding cells 
        ExCluster = [];
        for itWF = 1: length (WFs)
            if TP_latency(itWF) < 0.425 && awakeMeanRate(itWF) > 2 %&& max_burstIndex(itWF) < 0.1 %||  % this is letting some outliers in still add busrt Index cap? 
                InCluster1 = [InCluster1;itWF]; % narrow wf ins 
            elseif TP_latency(itWF) >= 0.425  && awakeMeanRate(itWF) > 1.5 && max_burstIndex(itWF) < 0.02 % 0.1 is what i used to make the final clusters 
                InCluster2 = [InCluster2;itWF]; %wide wf ins -          
            else
                ExCluster = [ExCluster;itWF]; 
            end            
        end 
        % making clusters for DG Excitatory and CA3 excitatory with the
        % same IN exclusion for now - maybe this needs to be changed later
        % since there might be different interneurons in the CA3
        % region(probbably not an issue since I didn't see many while
        % cutting) 
        DG_ExCluster = DG_cluster(ismember(DG_cluster, ExCluster));
        CA3_ExCluster = CA3_cluster(ismember(CA3_cluster, ExCluster)); %just save this variable 
        AMB_ExCluster = AMB_cluster(ismember(AMB_cluster, ExCluster)); 
 end

 function [wfPC1,wfPC2, pca_data_no_trail, diff_pca_data] = waveformPCA(DG_ExCluster,waveforms)
%makes wf-PCA components to be used in classificaiton
    %other options: 
    % 1) create a wave width at half height feature? - instead of the PCA on
    % the whole wave? - can get from extract DS2 artejact rejection 
    % 2) try just plotting sections of the second derrivative to see if you get
    %a bimodality 

    ex_waveforms = waveforms(DG_ExCluster);
    pca_data = [];
    for itEx = 1: length(ex_waveforms)
         if ~isnan(cell2mat(ex_waveforms(itEx))) %& length(cell2mat(ex_waveforms(itEx))) == 50  
             pca_data = [pca_data; interp1(1:50, cell2mat(ex_waveforms(itEx)),1:0.48:50,'spline')];
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
    
    pca_data_no_trail = padWF( : , firstPkInd:nSampsWF );%remove leading trail values because of zero padding 
    pca_data_no_trail = pca_data_no_trail(: , 6:end); %not sure why but there are zeros left in the first 5 columns and im taking them out like this. 
%     pca_data_no_trail = pca_data_aligned(:,maxShift+1:end);%data that
%     goes into the pca starts at the peak like S&B 
    
    % second derrivative 
    diff_pca_data = [];
    for itPCA = 1: size(pca_data_no_trail,1)
        diff_pca_data = [diff_pca_data; diff(pca_data_no_trail(itPCA,maxShift:end),2)]; %inserted 1:80 for the .8 ms time window % better time window - maxShift:end
    end
%     [coeff,score] =  pca(diff_pca_data(:,:)); %only hc number here checking full wf to see how it looks 
%     w_PC1 = coeff(:,1);
%     w_PC2 = coeff(:,2);
%     % weight of the PCA
%     wfPC1 = score(:,1);
%     wfPC2 = score(:,2);
    
    %load PCA coeffs from probe data 
    load ("/Users/isavars/Documents/phd_data/MatLabData/waveform_pca_coeffs_from_probes.mat", 'w_PC1', 'w_PC2')
    w_PC1 = w_PC1(1:size(diff_pca_data,2));
    w_PC2 = w_PC2(1:size(diff_pca_data,2));

    %lines to add instead of actually doing the pca 
    xCentred = diff_pca_data - nanmean(diff_pca_data,1);
    wfPC1 = xCentred/w_PC1';
    wfPC2 = xCentred/w_PC2';

    %k-means on the 1st two PCs of the wfPCA - not using for clustering 
%     
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
    log_normalized_awakeMeanRate = normalize(log1p(awakeMeanRate));

    % Apply log normalization to BI
    BI = data(:, 3);
    log_normalized_BI = normalize(log1p(BI));

    % Combine the log-normalized feature with the rest of the data
    transformed_data = [data(:, 1), log_normalized_awakeMeanRate, data(:, 3), data(:, 4)];% log_normalized_BI, data(:, 4)]; %, data(:, 5)];% ]; %

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

function DS2_orientations = get_DS2_orientations(spatData, DG_ExCluster,elePos, tets, option)
% are the spikes on the shank the cell was recorded on pointing up, mix of
% both or down. 
    
    %get tetrodes from channels for tetrode data only 
    chan_tets = zeros(64, 2);
    chan_tets(:, 1) = (1:64)';
    chan_tets(:, 2) = repelem(1:16, 4)';
    %get channels from elePos
    channels  = elePos.DS2_channels; 
    tets_from_chans = cell(height(elePos),1);
    for ii = 1: height(elePos)
        % Find the matching indices in the first column of chan_tets
        [~, indices] = ismember(channels(ii,:), chan_tets(:, 1));
        % Replace the values in tets_from_chans with the corresponding values from the second column of chan_tets
        indices = indices(indices ~=0);
        tets_from_chans{ii,:} = chan_tets(indices, 2)';
    end

    DS2_orientations = zeros(height(spatData),1);    
    DS2_orientations_option = []; %making these to look at population values on a histogram 

    %loop through spatData(DG_ExCluster) and find the DS2 amplitudes and make a distance from DS2 discrete measure per shank or section of one shank probe 
    
    for it_DG_Ex = DG_ExCluster'
        for it_ep = 1: height(elePos)
            if strcmp(spatData.dataset(it_DG_Ex),elePos.dataset(it_ep))
                [~ , tet] = ismember(tets(it_DG_Ex),tets_from_chans{it_ep,:}); %this needs have the position of the tetrode the cell was recored on in elepos column format            
                if tet ~=0 
                    if option == 1
                        DS2_median_amplitude_per_shank = median(elePos.DS2_max_amplitude(it_ep,tet)); 
                        % ^^ this finds the value per trial for tetrode the
                        % cell was recorded from - theres only one eeg per
                        % tetrode normally so median is a bit unecessary but
                        % leaving in case there are ever two eegs per tet 
                        if DS2_median_amplitude_per_shank > 0.4
                            DS2_orientation = 1;
                        elseif DS2_median_amplitude_per_shank <= 0.4 || DS2_median_amplitude_per_shank >= -0.2 
                            DS2_orientation = 2;
                        elseif DS2_median_amplitude_per_shank < -0.2
                            DS2_orientation = 3;
                        elseif isnan(DS2_median_amplitude_per_shank)
                            DS2_orientation = 4; 
                        end
                        DS2_orientations(it_DG_Ex) = DS2_orientation;
                        DS2_orientations_option = [DS2_orientations_option;DS2_median_amplitude_per_shank];
                    elseif option == 2 
                        %maybe plot first for whole dataset to come up with cuttoffs 
                        DS2_mean_amplitude_per_shank = median(elePos.DS2_peak_to_trough_amplitude(it_ep,tet));
                        if DS2_mean_amplitude_per_shank > 0.4
                            DS2_orientation = 1;
                        elseif DS2_mean_amplitude_per_shank <= 0.4 || DS2_mean_amplitude_per_shank >= -0.05 
                            DS2_orientation = 2;
                        elseif DS2_mean_amplitude_per_shank < -0.05
                            DS2_orientation = 3;
                        elseif isnan(DS2_mean_amplitude_per_shank)
                            DS2_orientation = 4; 
                        end
                        DS2_orientations(it_DG_Ex) = DS2_orientation;
                        DS2_orientations_option = [DS2_orientations_option;DS2_mean_amplitude_per_shank];
                    elseif option == 3 %fill in for slope - plot first to come up with cuttoffs 
                        DS2_mean_slope_per_shank = median(elePos.DS2_slope(it_ep,tet));
                        if DS2_mean_slope_per_shank < -50
                            DS2_orientation = 1;
                        elseif DS2_mean_slope_per_shank <= 0 || DS2_mean_slope_per_shank >= -50
                            DS2_orientation = 2;
                        elseif DS2_mean_slope_per_shank > 0
                            DS2_orientation = 3;
                        elseif isnan(DS2_mean_slope_per_shank)
                            DS2_orientation = nan; 
                        end
                        DS2_orientations(it_DG_Ex) = DS2_orientation;
                        DS2_orientations_option = [DS2_orientations_option;DS2_mean_slope_per_shank];
                    end 
                else 
                    DS2_orientations(it_DG_Ex) = nan;
                    DS2_orientations_option = [DS2_orientations_option;nan];
                end
            end
        end 
    end
    %histogram(DS2_orientations_option,'NumBins', 10);
end
% Function for jittering data points
function jitteredData = jitter(data)
    jitterRange = 0.5;  % Adjust as needed
    jitteredData = data + (rand(size(data))-0.5) * jitterRange;
end

% Function for jittering a categorical variable
function jitteredData = jitterCategorical(data)
    jitterRange = 0.3;  % Adjust as needed
    jitteredData = zeros(size(data));

    for i = 1:numel(data)
        jitter = (rand(1)-0.5) * jitterRange;
        if strcmp(data{i}, 'up')
            jitteredData(i) = 1 + jitter;
        elseif strcmp(data{i}, 'inverting')
            jitteredData(i) = 2 + jitter;
        elseif strcmp(data{i}, 'down')
            jitteredData(i) = 3 + jitter;
        else
            jitteredData(i) = 4 + jitter;
        end
    end
end

function [meanRate_per_tet,ratio_silent_or_active_per_tet]= make_silent_vs_active_per_tet(spatData, AMB_cluster,cellInfo, DG_ExCluster,amb_wakeMeanRate, sleep_idx, wake_idx )
%this needs to read in cells from amb cluster and provide say if they are active in sleep only or at least one wake env 
% then take the average mean firing rate of all cells on that tet and rank
% them from low to high - introdicting cuttoffs where the ratio of active
% in wake vs only active in sleep flips - that can be considered the hillus
% and then when it filps back it can be considered CA3 - go per rat and do
% all depths then see if it aligns with DS2 - can only spend time on this
% if you find a moment. 

    %itterate over cells and find if the cells are active in sleep only or
    %in run trials
    spatData_amb = spatData(DG_ExCluster,:);
    sleep_idx_amb = sleep_idx(DG_ExCluster);
    %spatData_amb = spatData(AMB_cluster,:);
    silent_or_active = zeros(height(spatData_amb),1);
    for itC = 1: height(spatData_amb)
        if all(spatData.nSpks(itC,wake_idx{itC}) < 75)
            silent_or_active(itC) = 0; %active in sleep only
        else 
            silent_or_active(itC) = 1; %active in wake trial
        end
    end
    %get the mean rate per tetrode 
    meanRate_per_tet = zeros(height(spatData_amb),1);
    ratio_silent_or_active_per_tet = zeros(height(spatData_amb),1);
    %cellInfo = cellInfo(AMB_cluster,:);
    cellInfo = cellInfo(DG_ExCluster,:);
    for jj=1:length(cellInfo)
        ID = cellInfo(jj, 1);
        tet = cellInfo(jj, 2);
        age = cellInfo(jj, 3);
        dates = cellInfo(jj, 5);      
        % Find indices of rows with the same ID, tet, age, and date
        matchingIndices = find(cellInfo(:,1) == ID & cellInfo(:,2) == tet & cellInfo(:,3) == age & cellInfo(:,5) == dates);
        %get mean rate per tet 
        meanRate_per_tet(jj) = nanmean(spatData_amb.meanRate(matchingIndices,sleep_idx_amb(jj)));
        ratio_silent_or_active_per_tet(jj) = sum(silent_or_active(matchingIndices) ==0)/sum(silent_or_active(matchingIndices) ==1); %(sum(silent_or_active(matchingIndices) ==0) + sum(silent_or_active(matchingIndices) ==1));%
        if ratio_silent_or_active_per_tet(jj) == inf
            ratio_silent_or_active_per_tet(jj) = 0;
        end
    end

    
    
    %rank the tetrodes by mean rate and compara against proportion of acive
    %vs silent - create DG vs CA3 labels for the tets 

    %makes plot - but need to refine 
    
    % Get unique tet values
    [uniqueTets, ~, tetIndices] = unique(meanRate_per_tet);
    
    % Initialize proportions arrays
    prop0 = zeros(size(uniqueTets));
    prop1 = zeros(size(uniqueTets));
    
    % Calculate proportions
    for i = 1:length(uniqueTets)
        currentTetIndices = tetIndices == i;
        totalForCurrentTet = sum(currentTetIndices);
        prop0(i) = sum(silent_or_active(currentTetIndices) == 0) / totalForCurrentTet;
        prop1(i) = sum(silent_or_active(currentTetIndices) == 1) / totalForCurrentTet;
    end
    
    % Create stacked bar plot with indices on x-axis
    bar(1:length(uniqueTets), [prop0, prop1], 1, 'stacked'); % Adjusted here for dimension alignment
    legend('silent', 'active');
    xlabel('Tet Index');
    ylabel('Proportion');
    title('Stacked Bar Plot of Proportions');
    
    % If you want to annotate the bars with the actual tet values
    for i = 1:length(uniqueTets)
        text(i, 0, num2str(uniqueTets(i)), 'Rotation', 90, 'HorizontalAlignment', 'right', 'FontSize', 8);
    end

    
end
function [normalized_data] = normalize_data_by_age(data, DG_ExCluster, spatData)
    %here we need to loop over age and make means of the exciatroy cell data for each age
    %bin subtract for each row then loop over data and create
    %normalized_data by subrtacting the mean for the right age bin from
    %each row 
    
    %get ages to make age bins
    cellInfo = getCellInfo(spatData);
    age = cellInfo(:,3);
    age = age(DG_ExCluster);

    postwean_data_idx = [];
    prewean_data_idx = [];
    adult_data_idx = []; %because its actually 3 age bins 
   for ii = 1: height(age)
        if age(ii) >= 21 && age(ii) <= 32
            postwean_data_idx = [postwean_data_idx; ii];
        elseif age(ii) >= 16
            prewean_data_idx = [prewean_data_idx; ii];
        elseif age(ii) == 40
            adult_data_idx = [adult_data_idx; ii];
        end
   end 

   %make means of age groups for normalizing the data 
    postwean_data_mean = mean(data(postwean_data_idx));
    prewean_data_mean = mean(data(prewean_data_idx));
    adult_data_mean = mean(data(adult_data_idx));
    
    normalized_data = data; % Initialize with original data

    % Normalize data by subtracting the mean of the right age bin
    normalized_data(postwean_data_idx,:) = data(postwean_data_idx,:) - postwean_data_mean;
    normalized_data(prewean_data_idx,:) = data(prewean_data_idx,:) - prewean_data_mean;
    normalized_data(adult_data_idx,:) = data(adult_data_idx,:) - adult_data_mean;

end