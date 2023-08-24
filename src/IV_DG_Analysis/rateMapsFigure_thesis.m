function rateMapsFigure_thesis (data, electrode_positions, clusters, writeDir)
%1. makes an array of axes to contain rate maps, autocorrelograms and waveforms of cells from my defined clusters
%2. organized from most spatial to least spatial 
%3. mantaining a number order (Labels)from the original table 
%4. Including waveform and AC from max channel

% TO DO: 
% 1. make more adaptable to different data sets when picking which rate map 
% to delete - i could just do this on illustrator. 
% 2. not adapted to tetrode data 


% obtaining variables from spatData Table (I want to replace this with a
% separate function called loadSpatData)

load (data, 'spatData')
load (electrode_positions, 'elePos')
load (clusters, 'PCA2_clusters', 'DG_ExCluster','CA3_ExCluster', 'InCluster1', 'InCluster2')

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env;
    rMap = spatData.rMap;
    SI_spat = spatData.SI_spat;
    cellID = spatData.cellID;
    SpkTs = spatData.SpkTs;
    waveforms = spatData.waveforms;
    nSpks = spatData.nSpks;
    trialDur = spatData.trialDur;
    

    meanRate(~(env=='fam' | env=='nov' | env=='diff' | env=='sleep')) = NaN; %waht is this for? 
    burstIndex(~(env=='fam' | env=='nov'| env=='diff'  | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='diff' | env=='sleep')) = NaN;

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

        
%prepare data for AC and waveforms  - this finds the one from the max wake
%trial i need to find the one from the max trial in general even if its
%sleep and for that i need an adaptable trial length feature. - create

    for itSp = 1: length (nSpks) 
        [~, maxSpksPos] = max(nSpks(itSp,wake_idx{itSp})); %change back to all
        STs(itSp,:) = SpkTs(itSp, maxSpksPos);
        WFs (itSp,:) = waveforms(itSp, maxSpksPos);
        trial_duration (itSp,:) = trialDur(itSp, maxSpksPos);
        max_meanRate(itSp,:) = meanRate (itSp, maxSpksPos);
        BIs (itSp,:) = burstIndex(itSp, maxSpksPos);
    end
    
%convert trial duraiton to integer value so it works in AC calculation
trial_duration = round(trial_duration);
    
%gather age data from cellInfo 
    cellInfo = getCellInfo(spatData);
    corecorded = cellInfo(:,4);   

%obtain clusters - load saved output from class_cells 

%     [DS2_amplitude,wfPC_clusters, DG_ExCluster, InCluster] = class_cells(data, electrode_positions);

%     AboveDS2 =[];
%     BelowDS2 =[];
%     for ii = 1: length(DS2_amplitude)
%         if ~isnan(DS2_amplitude(ii))
%             if DS2_amplitude(ii) <= 0 
%                 AboveDS2 = [AboveDS2;DG_ExCluster(ii)]; 
%             elseif DS2_amplitude(ii) > 0
%                 BelowDS2 = [BelowDS2;DG_ExCluster(ii)];     
%             end
%         end
%     end 
% 
%     clusters = {AboveDS2, BelowDS2};

    cluster1 =[];
    cluster2 =[];
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 1 %granule cells go here and will always be the first pages of maps 
                cluster1 = [cluster1;DG_ExCluster(ii)]; 
        elseif PCA2_clusters(ii) == 2
            cluster2 = [cluster2;DG_ExCluster(ii)];
        end
    end 

% %     cluster1 = low_narrow;
%     cluster1 = InCluster2;
%     cluster2 = InCluster2;
    
    clusters = {[139;185;186;187;207;210;219;229;233]};%, cluster2};%, cluster2};


% create ranking of spatiallity in cluster and arrange from most spatial to
% least spatial based on the SI_spat score. 

    for itC = 1:length(clusters)
        spatRank = nanmean(SI_spat(clusters{itC},wake_idx{itC}),2);
        SpatRankCluster = zeros(length(clusters{itC}),2); 
        SpatRankCluster(:,1) = clusters{itC};
        SpatRankCluster(:,2) = spatRank;
        SpatRankCluster = sortrows(SpatRankCluster,2, 'descend');
        clusters{itC} = SpatRankCluster(:,1);
    end
    

% makes 1 figure per cluster
    
    maxRowPerFig = 5;
    axRowCount = 1;


    %   makes labels 

    textContent = strcat((extractBefore (cellID, '_')),' ','P',(extractAfter (cellID, 'P')));

    fig_count = 0;
    
    for itC = 1:length(clusters)
        for it_clu = 1: length(clusters{itC}) 
            if it_clu == 1 || axRowCount > maxRowPerFig
                axRowCount = 1;
                hFig = gra_multiplot(maxRowPerFig, 7, 'figborder', [2 1 1 1]);
                axArr = getappdata(hFig, 'axesHandles' ); % makes the axes     
            end       
            for it_rm = 1: 5
                    gra_plotmap(rMap{clusters{itC}(it_clu),it_rm}, 'parent', axArr(axRowCount,it_rm)); % to put stuff in the axes. 
            end
%             spk_crosscorr(cell2mat(STs(clusters{itC}(it_clu))),'AC',0.001,0.3,900,'plot', axArr(axRowCount,6));% store these somewhere instead of making them 
            spk_crosscorr(cell2mat(STs(clusters{itC}(it_clu))),'AC',0.001,0.3,trial_duration(clusters{itC}(it_clu)),'plot', axArr(axRowCount,6));% store these somewhere instead of making them 
                axis(axArr(axRowCount,6),'tickaligned');
                axArr(axRowCount,6).XAxis.TickLabelRotation = 0;
                axArr(axRowCount,6).YAxis.TickLabelRotation = 0;
                %xlabel(axArr(axRowCount,6), 'Seconds');
                ylabel(axArr(axRowCount,6), num2str(BIs(clusters{itC}(it_clu))));

            plot(axArr(axRowCount,7), cell2mat(WFs(clusters{itC}(it_clu))));
                axis(axArr(axRowCount,7),[0 75 -100 200]); 
                axArr(axRowCount,7).XTickLabel = {};
                xlabel(axArr(axRowCount,7), '1.5 ms', 'FontSize', 12); %this needs to change when its tetrode data 
                hYLabel = ylabel(axArr(axRowCount,7), 'Voltage (μV)');
                set(hYLabel, 'Position', [90 1 1]);
            text (axArr(axRowCount,1),-50,23,textContent(clusters{itC}(it_clu)), 'FontSize', 12); 
            axRowCount = axRowCount + 1;

            % Save the figure once for every 5 rows
            if mod(it_clu, 5) == 0
                %save the figures with the right names
                if itC == 1
                    cluster = 'Putatve Granule';
                elseif itC == 2
                    cluster = 'Putative Mossy';
                end
                fig_count = fig_count +1;
%                 savefig(hFig, [writeDir '/' cluster ': Group ' num2str(fig_count) '.fig'])
            end
        end
        fig_count = 0; 
   end 

end


