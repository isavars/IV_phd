function rateMapsFigure ()
%1. will make an array of axes to contain rate maps from my defined clusters
%2. organize from most spatial to least spatial 3. mantaining a number order
%from the original table (try to print this) 4. find a way to remove extra trials
load ('r889_P18-P22_spatData_fin.mat', 'spatData')

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env;
    rMap = spatData.rMap;
    cellNo = spatData.cellNo;
    SI_spat = spatData.SI_spat;
    

    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
%     rMap(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;

    for it_gm = 1: length (meanRate)
        meanMeanRate(it_gm) = mean(meanRate(it_gm,:));
        meanBurstIndex(it_gm) = mean(burstIndex(it_gm,:));
    end

    cluster1 = [];
    cluster2 = [];
    cluster3 = [];
    cluster4 = [];
    
    for jj = 1: length(meanMeanRate)
        if  (meanMeanRate(jj) > 1) && (meanMeanRate(jj) < 10) && (meanBurstIndex(jj) < 0.05) %mystery cells
            cluster1(jj) = jj;
        elseif  (meanMeanRate(jj) > 0.01) && (meanMeanRate(jj) < 0.99)  && (meanBurstIndex(jj) > 0.02) && (meanBurstIndex(jj) < 0.15) %potential granule cells
            cluster2(jj) = jj;
        elseif (meanMeanRate(jj) > 0.02) && (meanBurstIndex(jj) > 0.15) %potential mossy cells
            cluster3(jj) = jj;
        else 
            cluster4(jj) = jj;
        end    
    end 
    cluster1 = nonzeros(cluster1);
    cluster2 = nonzeros(cluster2);
    cluster3 = nonzeros(cluster3);
    
    NewCluster1 = transpose(cluster1);
    NewCluster2 = transpose(cluster2);
    NewCluster3 = transpose(cluster3);

% create ranking of spatiallity in cluster and arrange from most spatial to
% least spatial based on the SI_spat score. 

    spatRank = [];
    for jj = NewCluster1
        spatRank(jj)= mean(SI_spat(jj,1:4));
    end
    spatRank = nonzeros(spatRank);
    SpatRankCluster1 = zeros(length(cluster1),2);
    SpatRankCluster1(:,1) = NewCluster1;
    SpatRankCluster1(:,2) = spatRank;
    SpatRankCluster1 = sortrows(SpatRankCluster1,2, 'descend');
    NewCluster1 = SpatRankCluster1(:,1);
    
    spatRank = [];
    for jj = NewCluster2
        spatRank(jj)= mean(SI_spat(jj,1:4));
    end
    spatRank = nonzeros(spatRank);
    SpatRankCluster = zeros(length(cluster2),2);
    SpatRankCluster(:,1) = NewCluster2;
    SpatRankCluster(:,2) = spatRank;
    SpatRankCluster = sortrows(SpatRankCluster,2, 'descend');
    NewCluster2 = SpatRankCluster(:,1);
    
    spatRank = [];
    for jj = NewCluster3
        spatRank(jj)= mean(SI_spat(jj,1:4));
    end
    spatRank = nonzeros(spatRank);
    SpatRankCluster = zeros(length(cluster3),2);
    SpatRankCluster(:,1) = NewCluster3;
    SpatRankCluster(:,2) = spatRank;
    SpatRankCluster = sortrows(SpatRankCluster,2, 'descend');
    NewCluster3 = SpatRankCluster(:,1);
    

% 1 figure per cluster

% %Mystery cells    
%     hFig = gra_multiplot( length(cluster1), 4, 'figborder', [2 1 1 1] );  % will have to personalize fig borders to different environment 
%     axArr = getappdata(hFig, 'axesHandles' ); % makes the axes
% 
%     for it_clu = 1: length(cluster1) 
%         for it_rm = 1: 4
%             gra_plotmap(rMap{NewCluster1(it_clu),it_rm}, 'parent', axArr(it_clu,it_rm) ); % to put stuff in the axes. 
%         end
%     end     
% %possible granule cells    
%     hFig = gra_multiplot( length(cluster2), 4, 'figborder', [2 1 1 1] );  % will have to personalize fig borders to different environment 
%     axArr = getappdata(hFig, 'axesHandles' ); % makes the axes
% 
%     for it_clu = 1: length(cluster2) 
%         for it_rm = 1: 4
%             gra_plotmap(rMap{NewCluster2(it_clu),it_rm}, 'parent', axArr(it_clu,it_rm) ); % to put stuff in the axes. 
%         end
%     end
%possible mossy cells    
    hFig = gra_multiplot( length(cluster3), 4, 'figborder', [0.5 0.5 0.5 0.5] );  % will have to personalize fig borders to different environment 
    axArr = getappdata(hFig, 'axesHandles' ); % makes the axes

    for it_clu = 1: length(cluster3) 
        for it_rm = 1: 4
            gra_plotmap(rMap{NewCluster3(it_clu),it_rm}, 'parent', axArr(it_clu,it_rm) ); % to put stuff in the axes. 
        end
    end  
end