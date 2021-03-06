function rateMapsFigure ()
%1. makes an array of axes to contain rate maps, autocorrelograms and waveforms of cells from my defined clusters
%2. organized from most spatial to least spatial 
%3. mantaining a number order (Labels)from the original table 
%4. Including waveform and AC from max channel
% 5. make more adaptable to different data sets had to change 1:3 to 2:4 -
% add indecexes . Also when picking which rate map to delete - i
% could just do this on illustrator. 

load ('all_Adult_spatData_CDB.mat', 'spatData')

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
    

    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    
%make indexes for environment to use in comparisons with any file 
        maxFam = 2;
        FamInd = nan(length(meanRate),maxFam); 
        FamIndT = [];     
        NovInd = [];
        famCount = 0;
        
        for itCell= 1: length(meanRate)
            for itTrial = 1: 5
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
    
%make rate change score
    
    for it_gma = 1: length (meanRate)
        awakeMeanRate(it_gma,:) = mean(meanRate(it_gma,(FamInd(1):FamInd(2))));
    end

    sleepMeanRate = meanRate (:,5);
    
    rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));
    
%gather age data from cellInfo 

    cellInfo = getCellInfo();
    corecorded = cellInfo(:,4);   
%obtain clusters 
    
    for it_gm = 1: length (meanRate)
        meanMeanRate(it_gm) = mean(meanRate(it_gm,5));
        meanBurstIndex(it_gm) = mean(burstIndex(it_gm,5));
    end

    cluster1 = [];
    cluster2 = [];
    cluster3 = [];
    cluster4 = [];

        for jj = 1: length(meanMeanRate)
            if  (meanMeanRate(jj) >= 10)   %interneurons (meanMeanRate(jj) >= 1) && (meanMeanRate(jj) <= 10) && (meanBurstIndex(jj) <= 0.05)
                cluster1 = [cluster1; jj];
            elseif  (rateChange(jj) <= 0.3) && (meanMeanRate(jj) >= 0.01) && (meanMeanRate(jj) <= 1.5)  && (meanBurstIndex(jj) >= 0.02) && (meanBurstIndex(jj) < 0.15) %potential granule cells
                cluster2 = [cluster2; jj];
            elseif (meanMeanRate(jj) <= 1.5) && (meanBurstIndex(jj) >= 0.15) && (corecorded(jj) > 3) %potential CA3 cells
                cluster4 = [cluster4; jj];
            elseif (meanMeanRate(jj) >= 0.02) && (meanBurstIndex(jj) >= 0.15) %potential mossy cells
                cluster3 = [cluster3; jj];

%             if  (meanMeanRate(jj) > 1) && (meanMeanRate(jj) < 10) && (meanBurstIndex(jj) < 0.05) %mystery cells
%                 cluster1 = [cluster1; jj];
%             elseif  (rateChange(jj) <= 0.3) && (meanMeanRate(jj) > 0.01) && (meanMeanRate(jj) < 1.5)  && (meanBurstIndex(jj) > 0.02) && (meanBurstIndex(jj) < 0.15) %potential granule cells
%                 cluster2 = [cluster2; jj];
%             elseif (meanMeanRate(jj) > 0.02) && (meanBurstIndex(jj) > 0.15) %potential mossy cells
%                 cluster3 = [cluster3; jj];
%             else 
%                 cluster4 = [cluster4; jj];
            end    
        end 
    
    NewCluster1 = transpose(cluster1);
    NewCluster2 = transpose(cluster2);
    NewCluster3 = transpose(cluster3);
    NewCluster4 = transpose(cluster4);
    

% create ranking of spatiallity in cluster and arrange from most spatial to
% least spatial based on the SI_spat score. 

    spatRank = [];
    for jj = NewCluster1
        spatRank(jj)= mean(SI_spat(jj,1:3));
    end
    spatRank = nonzeros(spatRank);
    SpatRankCluster1 = zeros(length(cluster1),3); %made changes here 13/04/2020
    SpatRankCluster1(:,1) = NewCluster1;
    SpatRankCluster1(:,2) = spatRank;
    SpatRankCluster1 = sortrows(SpatRankCluster1,2, 'descend');
    NewCluster1 = SpatRankCluster1(:,1);
    
    spatRank = [];
    for jj = NewCluster2
        spatRank(jj)= mean(SI_spat(jj,1:3));
    end
    spatRank = nonzeros(spatRank);
    SpatRankCluster = zeros(length(cluster2),2);
    SpatRankCluster(:,1) = NewCluster2;
    SpatRankCluster(:,2) = spatRank;
    SpatRankCluster = sortrows(SpatRankCluster,2, 'descend');
    NewCluster2 = SpatRankCluster(:,1);
    
    spatRank = [];
    for jj = NewCluster3
        spatRank(jj)= mean(SI_spat(jj,1:3));
    end
    spatRank = nonzeros(spatRank);
    SpatRankCluster = zeros(length(cluster3),2);
    SpatRankCluster(:,1) = NewCluster3;
    SpatRankCluster(:,2) = spatRank;
    SpatRankCluster = sortrows(SpatRankCluster,2, 'descend');
    NewCluster3 = SpatRankCluster(:,1);
    
    spatRank = [];
    for jj = NewCluster4
        spatRank(jj)= mean(SI_spat(jj,1:3));
    end
    spatRank = nonzeros(spatRank);
    SpatRankCluster = zeros(length(cluster4),2);
    SpatRankCluster(:,1) = NewCluster4;
    SpatRankCluster(:,2) = spatRank;
    SpatRankCluster = sortrows(SpatRankCluster,2, 'descend');
    NewCluster4 = SpatRankCluster(:,1);
    

% makes 1 figure per cluster
    
    maxRowPerFig = 5;
    axRowCount = 1;
    
    %   makes labels 

    textContent = cellID; %IV 04/12/21 strcat((extractBefore (cellID, '_')),'P',(extractAfter (cellID, 'P')));
    
%     for it_clu = 1: length(NewCluster1) 
%         if it_clu == 1 || axRowCount > maxRowPerFig
%             axRowCount = 1;
%             hFig = gra_multiplot(maxRowPerFig, 6, 'figborder', [2 1 1 1]);
%             axArr = getappdata(hFig, 'axesHandles' ); % makes the axes     
%         end       
%         for it_rm = 1: 4
%                 gra_plotmap(rMap{NewCluster1(it_clu),it_rm}, 'parent', axArr(axRowCount,it_rm)); % to put stuff in the axes. 
% %                 if env(NewCluster1(it_clu),1) ~= env(NewCluster1(it_clu),4) 
% %                     map= gra_plotmap(rMap{NewCluster1(it_clu),4}, 'parent', axArr(axRowCount, 4),'text_pos', 'none');
% %                     delete(map)
% %                 end
%         end
% %         spk_crosscorr(cell2mat(STs(NewCluster1(it_clu))),'AC',0.002,0.1,900,'plot', axArr(axRowCount,5));% store these somewhere instead of making them 
%         plot(axArr(axRowCount,6), cell2mat(WFs(NewCluster1(it_clu))));
%             axis(axArr(axRowCount,6),[0 50 -100 100]); 
%         text (axArr(axRowCount,1),-50,23,textContent(NewCluster1(it_clu)), 'FontSize', 16); 
%         axRowCount = axRowCount + 1;
%     end 

   for it_clu = 1: length(NewCluster2) 
        if it_clu == 1 || axRowCount > maxRowPerFig
            axRowCount = 1;
            hFig = gra_multiplot(maxRowPerFig, 6, 'figborder', [2 1 1 1] );
            axArr = getappdata(hFig, 'axesHandles' ); % makes the axes 
        end
        for it_rm = 1: 4
            gra_plotmap(rMap{NewCluster2(it_clu),it_rm}, 'parent', axArr(axRowCount,it_rm) ); % to put stuff in the axes. 
%              if env(NewCluster2(it_clu),1) ~= env(NewCluster2(it_clu),4) 
%                     map= gra_plotmap(rMap{NewCluster2(it_clu),4}, 'parent', axArr(axRowCount, 4),'text_pos', 'none');
%                     delete(map)
%              end
        end
        spk_crosscorr(cell2mat(STs(NewCluster2(it_clu))),'AC',0.002,0.3,900,'plot', axArr(axRowCount,5));% store these somewhere instead of making them 
        plot(axArr(axRowCount,6), cell2mat(WFs(NewCluster2(it_clu))));
            axis(axArr(axRowCount,6),[0 50 -100 100]);
        text (axArr(axRowCount,1),-50,23,textContent(NewCluster2(it_clu)), 'FontSize', 16); 
        axRowCount = axRowCount + 1;
    end 

%     for it_clu = 1: length(NewCluster3) 
%         if it_clu == 1 || axRowCount > maxRowPerFig
%             axRowCount = 1;
%             hFig = gra_multiplot(maxRowPerFig, 6, 'figborder', [2 1 1 1] );
%             axArr = getappdata(hFig, 'axesHandles' ); % makes the axes 
%         end
%         for it_rm = 1: 4
%             gra_plotmap(rMap{NewCluster3(it_clu),it_rm}, 'parent', axArr(axRowCount,it_rm) ); % to put stuff in the axes. 
% %              if env(NewCluster3(it_clu),1) ~= env(NewCluster3(it_clu),4) 
% %                     map= gra_plotmap(rMap{NewCluster3(it_clu),4}, 'parent', axArr(axRowCount, 4),'text_pos', 'none');
% %                     delete(map)
% %              end            
%         end
%          spk_crosscorr(cell2mat(STs(NewCluster3(it_clu))),'AC',0.002,0.3,900,'plot', axArr(axRowCount,5));% store these somewhere instead of making them 
%         plot(axArr(axRowCount,6), cell2mat(WFs(NewCluster3(it_clu))));
%             axis(axArr(axRowCount,6),[0 50 -100 100]);
%         text (axArr(axRowCount,1),-50,23,textContent(NewCluster3(it_clu)), 'FontSize', 16); 
%         axRowCount = axRowCount + 1;
%     end 
    
%     for it_clu = 1: length(NewCluster4) 
%         if it_clu == 1 || axRowCount > maxRowPerFig
%             axRowCount = 1;
%             hFig = gra_multiplot(maxRowPerFig, 6, 'figborder', [2 1 1 1] );
%             axArr = getappdata(hFig, 'axesHandles' ); % makes the axes 
%         end
%         for it_rm = 1: 4
%             gra_plotmap(rMap{NewCluster4(it_clu),it_rm}, 'parent', axArr(axRowCount,it_rm) ); % to put stuff in the axes. 
% %              if env(NewCluster4(it_clu),1) ~= env(NewCluster4(it_clu),4) 
% %                     map= gra_plotmap(rMap{NewCluster4(it_clu),4}, 'parent', axArr(axRowCount, 4),'text_pos', 'none');
% %                     delete(map)
% %              end  
% 
%         end
%          spk_crosscorr(cell2mat(STs(NewCluster4(it_clu))),'AC',0.002,0.3,900,'plot', axArr(axRowCount,5));% store these somewhere instead of making them 
%         plot(axArr(axRowCount,6), cell2mat(WFs(NewCluster4(it_clu))));
%             axis(axArr(axRowCount,6),[0 50 -100 100]);
%         text (axArr(axRowCount,1),-50,23,textContent(NewCluster4(it_clu)), 'FontSize', 16); 
%         axRowCount = axRowCount + 1;
%     end 
 
end