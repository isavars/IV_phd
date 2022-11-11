function rateMapsFigure_SfN (spatData)
%1. makes an array of axes to contain rate maps, autocorrelograms and waveforms of cells from my defined clusters
%2. organized from most spatial to least spatial 
%3. mantaining a number order (Labels)from the original table 
%4. Including waveform and AC from max channel
% 5. make more adaptable to different data sets had to change 1:3 to 2:4 -
% add indecexes . Also when picking which rate map to delete - i
% could just do this on illustrator. 

%load ('spatData_r1099.mat', 'spatData')

% obtaining variables from spatData Table (I want to replace this with a
% separate function called loadSpatData)

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env;
    rMap = spatData.rMap;
    %cellNo = spatData.cellNo;
    SI_spat = spatData.SI_spat;
    cellID = spatData.cellID;
    SpkTs = spatData.SpkTs;
    waveforms = spatData.waveforms;
    nSpks = spatData.nSpks;
    

    meanRate(~(env=='fam' | env=='nov' | env=='diff' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov'| env=='diff'  | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='diff' | env=='sleep')) = NaN;
    
%make indexes for environment to use in comparisons with any file 
        maxFam = 2;
        FamInd = nan(length(meanRate),maxFam); 
        FamIndT = [];     
        NovInd = [];
        famCount = 0;
        
        for itCell= 1: length(meanRate)
            for itTrial = 1: 6
                if contains(cast(spatData.env(itCell,itTrial),'char'),'fam')
                    FamIndT(itCell,itTrial) = itTrial;
                    famCount = famCount + 1;
                    if famCount <= maxFam
                         FamInd(itCell,famCount)= transpose(nonzeros(FamIndT(itCell,itTrial)));
                    end
                elseif strcmp(cast(spatData.env(itCell,itTrial),'char'),'nov')
                    NovInd(itCell,itTrial) = itTrial;
                    NovInd = nonzeros(NovInd);
                end 
            end
            famCount = 0;
        end
        
%prepare data for AC and waveforms  

    for itSp = 1: length (nSpks) 
        [~, maxSpksPos] = max(nSpks(itSp,2:4));
        STs(itSp,:) = SpkTs(itSp, maxSpksPos);
        WFs (itSp,:) = waveforms(itSp, maxSpksPos);
    end
    
    
%gather age data from cellInfo 

    cellInfo = getCellInfo(spatData);
    corecorded = cellInfo(:,4);   
%obtain clusters 

    [cluster1,cluster2,cluster3,cluster4] = makeClusters(spatData);
    
    NewCluster1 = transpose(cluster1);
    NewCluster2 = transpose(cluster2);
    NewCluster3 = transpose(cluster3);
    NewCluster4 = transpose(cluster4);
    
    clusters = {NewCluster1,NewCluster2,NewCluster3,NewCluster4}; %array of cluster names  

% create ranking of spatiallity in cluster and arrange from most spatial to
% least spatial based on the SI_spat score. 

    for itC = 1:length(clusters)
        spatRank = [];
        for jj = clusters{itC}
            spatRank(jj)= mean(SI_spat(jj,1:5));
        end
        spatRank = nonzeros(spatRank);
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
    
    for itC = 1:length(clusters)
        for it_clu = 1: length(clusters{itC}) 
            if it_clu == 1 || axRowCount > maxRowPerFig
                axRowCount = 1;
                hFig = gra_multiplot(maxRowPerFig, 7, 'figborder', [2 1 1 1]);
                axArr = getappdata(hFig, 'axesHandles' ); % makes the axes     
            end       
            for it_rm = 1: 5
                    gra_plotmap(rMap{clusters{itC}(it_clu),it_rm}, 'parent', axArr(axRowCount,it_rm)); % to put stuff in the axes. 
    %                 if env(NewCluster1(it_clu),1) ~= env(NewCluster1(it_clu),4) 
    %                     map= gra_plotmap(rMap{NewCluster1(it_clu),4}, 'parent', axArr(axRowCount, 4),'text_pos', 'none');
    %                     delete(map)
    %                 end
            end
            spk_crosscorr(cell2mat(STs(clusters{itC}(it_clu))),'AC',0.002,0.3,900,'plot', axArr(axRowCount,6));% store these somewhere instead of making them 
            plot(axArr(axRowCount,7), cell2mat(WFs(clusters{itC}(it_clu))));
                axis(axArr(axRowCount,7),[0 50 -100 100]); 
            text (axArr(axRowCount,1),-50,23,textContent(clusters{itC}(it_clu)), 'FontSize', 16); 
            axRowCount = axRowCount + 1;
        end
   end 


end