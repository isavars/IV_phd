function rateMapsFigure_thesis (data, electrode_positions, clusters, writeDir)
%1. makes an array of axes to contain rate maps, autocorrelograms and waveforms of cells from my defined clusters
%2. organized from most spatial to least spatial 
%3. mantaining a number order (Labels)from the original table 
%4. Including waveform and AC from max channel

% TO DO: 
% 1. make more adaptable to different data sets when picking which rate map 
% to delete - i could just do this on illustrator. 
% 2. add a trial duration feature (this will be in spat data) to use
% in the AC making bit when its not 900s long (so for sleep) 


% obtaining variables from spatData Table (I want to replace this with a
% separate function called loadSpatData)

load (data, 'spatData')
load (electrode_positions, 'elePos')
load (clusters, 'PCA2_clusters', 'DG_ExCluster')

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env;
    rMap = spatData.rMap;
    SI_spat = spatData.SI_spat;
    cellID = spatData.cellID;
    SpkTs = spatData.SpkTs;
    waveforms = spatData.waveforms;
    nSpks = spatData.nSpks;
    

    meanRate(~(env=='fam' | env=='nov' | env=='diff' | env=='sleep')) = NaN; %waht is this for? 
    burstIndex(~(env=='fam' | env=='nov'| env=='diff'  | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='diff' | env=='sleep')) = NaN;

        
%prepare data for AC and waveforms  - this finds the one from the max wake
%trial i need to find the one from the max trial in general even if its
%sleep and for that i need an adaptable trial length feature. - create

    for itSp = 1: length (nSpks) 
        [~, maxSpksPos] = max(nSpks(itSp,1:5));
        STs(itSp,:) = SpkTs(itSp, maxSpksPos);
        WFs (itSp,:) = waveforms(itSp, maxSpksPos);
        max_meanRate(itSp,:) = meanRate (itSp, maxSpksPos);
    end
    
    
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
%     
%             end
%         end
%     end 
% 
%     clusters = {AboveDS2, BelowDS2};

    cluster1 =[];
    cluster2 =[];
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 2 %granule cells go here and will always be the first pages of maps 
%             if max_meanRate(ii) < 0.5 %temporary rate filter for granule cell cluster 
                cluster1 = [cluster1;DG_ExCluster(ii)]; 
%             else 
%             end
        elseif PCA2_clusters(ii) == 1
            cluster2 = [cluster2;DG_ExCluster(ii)];

        end
    end 
    
    clusters = {cluster1, cluster2};


% create ranking of spatiallity in cluster and arrange from most spatial to
% least spatial based on the SI_spat score. 

    for itC = 1:length(clusters)
        spatRank = nanmean(SI_spat(clusters{itC},1:5),2);
        SpatRankCluster = zeros(length(clusters{itC}),2); 
        SpatRankCluster(:,1) = clusters{itC};
        SpatRankCluster(:,2) = spatRank;
        SpatRankCluster = sortrows(SpatRankCluster,2, 'descend');
        clusters{itC} = SpatRankCluster(:,1);
    end
    

% makes 1 figure per cluster
    
    maxRowPerFig = 6;
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
            spk_crosscorr(cell2mat(STs(clusters{itC}(it_clu))),'AC',0.001,0.3,900,'plot', axArr(axRowCount,6));% store these somewhere instead of making them 
                axis(axArr(axRowCount,7),[0 75 -100 200]);
            plot(axArr(axRowCount,7), cell2mat(WFs(clusters{itC}(it_clu))));
                axis(axArr(axRowCount,7),[0 75 -100 200]); 
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
                savefig(hFig, [writeDir '/' cluster ': Group ' num2str(fig_count) '.fig'])
            end
        end
        fig_count = 0; 
   end 


end