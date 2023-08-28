function get_DS2_and_theta_wrt_SpkTs (data, electrodes, clusters)
% I want this function to do two things 1. compare cell spike times to DS2
% spike times and create a peri-stimulus time histogram. this requires
% spk_psth (loop over cells and make the overal mean) - then do stats on
% this. 2. compare spike times to peak theta frequencey (get peak frequency 
% with eeg_powerspec) to find out if they have a phase preference uses 
% eeg_instfrequency. plot this somehow and do
% stats. 

%big problem - the spike times for ds2 are based on shortened trials and
%didnt preserve the oridingal trial length when they were cut down to the
%sleep epochs - youb are going to have to reverse engineer this - the sleep
%epochs are found in all_sleep_data_for_DS2_sleepData.mat basically youll
%need to locate which epoch each DS2 spike time is in and add the start
%time of the epoch to that timestamp. 

    % load spike times for each cluster during the sleep trial 
    load(data, 'spatData')

    %get sleep spike times    
    SpkTs_sleep = cell(height(spatData),1);
    for ii = 1:height(spatData)
        sleep_idx = find(strcmp(string(spatData.env(ii,:)), 'sleep'));
        SpkTs_sleep{ii} = spatData.SpkTs{ii,sleep_idx};
    end
    
    % all these spike times need to be adjusted back to cut trial times
    % using sleep epochs :( 
    %load sleep epochs
    %msubfunciton that shifts the time stamps using sleep epochs -
    %only needs to be done for DS2 comparison  
    [adjusted_SpkTs_sleep] = adjust_timestamps(spatData,SpkTs_sleep);
    %get DS2 spiketimes per cell 
    [DS2_SpkTs] =get_DS2_Spks(spatData, electrodes);

    
    %load clusters and get spike times for each cluster     
    %load clusters
    load(clusters, 'PCA2_clusters', 'DG_ExCluster','CA3_ExCluster')     
    %make clusters from PCA2_clusters
    mossy_cluster =[];
    granule_cluster =[]; 
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 2
            mossy_cluster = [mossy_cluster;DG_ExCluster(ii)]; 
        elseif PCA2_clusters(ii) == 1
            granule_cluster = [granule_cluster;DG_ExCluster(ii)];
        end
    end 

    pyramidal_cluster = CA3_ExCluster;

%     granule_cluster = [4; 8; 10; 11; 40; 86; 87; 235; 240]; %spatial only


    clusters = {granule_cluster, mossy_cluster, pyramidal_cluster}; 
    clusternames = {'granule', 'mossy', 'CA3'};


    %loop throught all the cells and get the peri-stimulus time histogram
    %for each. Probably make this its own function 


    %make peri-stimulus time histogram - trying for 1st cell
    spikeTrain = adjusted_SpkTs_sleep; %spike times in seconds
    DS2Times= DS2_SpkTs;
    window = [-250 250]; %1s in ms
    binSize = 20;


    %loop over each cluster and make means for each cell time 
    for it_clu = 1:length(clusters)
        cluster = clusters{it_clu};
        clustername = clusternames{it_clu};
        allCellsMeanHists = []; % Initialize
        allCellsStdHists = [];  % Initialize
        allCellsPeakRatios = []; % Initialize
    
        for i = cluster'
            [meanHist, stdHist, binEdgeTimes] = spk_psth(spikeTrain{i}, DS2Times{i}, window, binSize); % Assume you modified your function to also return stdHist
    
            allCellsMeanHists = [allCellsMeanHists; meanHist];
            allCellsStdHists = [allCellsStdHists; stdHist];

           
            % Define baseline window indices (for example, considering the first half)
            baselineIndices = 1:floor(length(binEdgeTimes)/2);
            
            % Calculate baseline rate
            baselineRate = nanmean(meanHist(baselineIndices));

            peakRate
            
            % Calculate the firing rate ratio for each bin
            firingRateRatios = meanHist / baselineRate;
    
            allCellsPeakRatios = [allCellsPeakRatios; firingRateRatios];
        end
        
        populationMean = nanmean(allCellsMeanHists, 1);
        populationStd = nanstd(allCellsMeanHists, [], 1);
        
        figure;
        bar(binEdgeTimes, populationMean, 'histc');
        title(['Population Mean Response:' clustername] );
        hold on;
%         fill([binEdgeTimes, fliplr(binEdgeTimes)], [populationMean+populationStd, fliplr(populationMean-populationStd)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        figure;
        plot(binEdgeTimes, populationMean, 'LineWidth', 2);
        xlabel('Time (s)');
        ylabel('Mean Spike Rate');
        title(['Population Mean Response:' clustername] );
        xlim([min(binEdgeTimes) max(binEdgeTimes)]);
        ax = gca;
        ax.XTick = [min(binEdgeTimes) 0 max(binEdgeTimes)]; % Set x-axis ticks to show start, middle (stimulus time), and end of the window
        grid on;

        populationPeakRatioMean = nanmean(allCellsPeakRatios, 1);

        % Plotting peak firing ratio
        figure;
        plot(binEdgeTimes, populationPeakRatioMean, 'LineWidth', 2);
        xlabel('Time (s)');
        ylabel('Mean Peak Firing Ratio');
        title(['Mean Peak Firing Ratio for ' clustername]);
        grid on;


        
    end

end
    
function [DS2_SpkTs] =get_DS2_Spks(spatData, electrodes)

    %load elepos
    load(electrodes, 'elePos')
    %make_DS2_SpkTs
    DS2_SpkTs = cell(height(spatData), 1);
    
    for it_C = 1: height(spatData)
        currentDataset = spatData.dataset(it_C);
        it_ep = find(strcmp(currentDataset, elePos.dataset), 1, 'first'); % find the matching index        
        if isempty(it_ep)
            DS2_SpkTs{it_C} = [];
            continue;
        end        
        if ~isempty(elePos.DS2_spike_times{it_ep}) && size(elePos.DS2_spike_times{it_ep}, 1) > 0
            %open each cell and time-adjust the spike times - save in
            %seconds
            DS2_SpkTs{it_C} = elePos.DS2_spike_times{it_ep}{1, 1};
            temp_SpkTs = cell2mat(DS2_SpkTs{it_C})./10^6;
            DS2_SpkTs{it_C} = temp_SpkTs;
        else
            DS2_SpkTs{it_C} = [];
        end
    end

end     

function [adjusted_SpkTs_sleep] = adjust_timestamps(spatData,SpkTs_sleep);
    load('all_sleep_data_for_DS2_sleepData.mat', 'sleepData')
    adjusted_SpkTs_sleep = cell(size(SpkTs_sleep)); % Initialize with the same size as SpkTs_sleep
    
    % Process each dataset's spike times
    for it_C = 1: length(SpkTs_sleep)
        currentDataset = spatData.dataset(it_C);
        sleepData_idx = find(strcmp(currentDataset, sleepData.dataset), 1, 'first'); % Find the matching index
        current_SpkT = SpkTs_sleep{it_C};
    
        % Get the SWS_epochs for the current dataset
        current_SWS_epochs = sleepData.SWS_epochs{sleepData_idx}; % Assuming that SWS_epochs is a cell array where each cell contains the 2xN double for a dataset
    
        starts = current_SWS_epochs(1, :);
        stops = current_SWS_epochs(2, :);
        adjusted_SpkT = current_SpkT; % Copy of current spike times to adjust
    
        for it_SpkTs = 1: length(current_SpkT)
            spikeTime = current_SpkT(it_SpkTs);
            cumulativeTime = 0; % To keep track of the time in previous epochs
    
            % Iterate over each epoch to check if spikeTime is within any epoch
            for epoch_idx = 1:length(starts)
                startTime = starts(epoch_idx);
                stopTime = stops(epoch_idx);
    
                % If spike is within this epoch, adjust and break
                if spikeTime >= startTime && spikeTime <= stopTime
                    adjusted_SpkT(it_SpkTs) = spikeTime - startTime + cumulativeTime;
                    break; 
                else
                    % If not, add this epoch's duration to the cumulative time
                    cumulativeTime = cumulativeTime + (stopTime - startTime);
                end
            end
        end
        
        % Store adjusted spike timestamps in adjusted_SpkTs_sleep
        adjusted_SpkTs_sleep{it_C} = adjusted_SpkT;
    end
    
end       
            