function BurstIndexVsMeanRatePlot ()
    
load ('r889_P18-P22_spatData_fin.mat', 'spatData')


    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
%     corecorded = activeCells();
    corecorded = corecordedCells(spatData);
    env = spatData.env; 
   
    
    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
   

%one scatterplot per trial pooled across days 
%     for jj= 1:5
%         figure()
%         scatter(meanRate(:,jj),burstIndex(:,jj));
%             plotTitle = strcat({'Burst Index Vs Mean Rate'}, {':'}, string(env(1,jj)));
%             title(plotTitle)
%             xlabel('Mean Rate (seconds)')
%             set(gca, 'xscale', 'log')
%             ylabel('Burst Index (%)')
%     end 
%     

    cluster1 = [];
    cluster2 = [];
    cluster3 = [];
    cluster4 = [];
    
    for jj = 1: length(meanRate(:,1))
        if  (meanRate(jj) > 1) && (meanRate(jj) < 10) && (burstIndex(jj) < 0.05) %mystery cells
            cluster1(jj) = jj;
        elseif (meanRate(jj) < 0.8) && (meanRate(jj) > 0.01) && (burstIndex(jj) < 0.15) && (burstIndex(jj) > 0.02) %potential granule cells
            cluster2(jj) = jj;
        elseif (meanRate(jj) > 0.02) && (burstIndex(jj) > 0.15) %potential mossy cells
            cluster3(jj) = jj;
        else 
            cluster4(jj) = jj;
        end    
    end 
    cluster1 = nonzeros(cluster1);
    cluster2 = nonzeros(cluster2);
    cluster3 = nonzeros(cluster3);
    
 

 % one scatterplot for all trials pooled across days 
 
%     meanRate = meanRate (:);
%     burstIndex = burstIndex (:);
%     
%     figure()
%     scatter(meanRate(:),burstIndex(:))
%         title('Burst Index Vs Mean Rate: All trials')
%         xlabel('Mean Rate (seconds)')
%         set(gca, 'xscale', 'log')
%         ylabel('Burst Index')

% scatterplots for the mean and max mean rate and burst index (all trials)
    
    for it_gm = 1: length (meanRate)
        meanMeanRate(it_gm) = mean(meanRate(it_gm,:));
        meanBurstIndex(it_gm) = mean(burstIndex(it_gm,:));
    end

    for it_gx = 1: length (meanRate)
        maxMeanRate(it_gx) = max(meanRate(it_gx,:));
        maxBurstIndex(it_gx) = max(burstIndex(it_gx,:));
    end    
    
    figure()
    scatter(meanMeanRate,meanBurstIndex)
        title('Means of Burst Index Vs Mean Rate: All trials')
        xlabel('Mean Rate (seconds)')
        set(gca, 'xscale', 'log')
        ylabel('Burst Index')

    figure()
    scatter(maxMeanRate,maxBurstIndex)
        title('Max values of Burst Index Vs Mean Rate: All trials')
        xlabel('Mean Rate (seconds)')
        set(gca, 'xscale', 'log')
        ylabel('Burst Index')

% one scatterplot for the max and mean rate and burst index for awake only and
% sleep only 

    for it_gma = 1: length (meanRate)
        meanMeanRate(it_gma) = mean(meanRate(it_gma,1:4));
        meanBurstIndex(it_gma) = mean(burstIndex(it_gma,1:4));
    end

    for it_gxa = 1: length (meanRate)
        maxMeanRate(it_gxa) = max(meanRate(it_gxa,1:4));
        maxBurstIndex(it_gxa) = max(burstIndex(it_gxa,1:4));
    end    
    
    figure()
    scatter(meanMeanRate,meanBurstIndex)
        title('Means of Burst Index Vs Mean Rate: Awake Trials')
        xlabel('Mean Rate (seconds)')
        set(gca, 'xscale', 'log')
        ylabel('Burst Index')

    figure()
    scatter(maxMeanRate,maxBurstIndex)
        title('Max values of Burst Index Vs Mean Rate: Awake Trials')
        xlabel('Mean Rate (seconds)')
        set(gca, 'xscale', 'log')
        ylabel('Burst Index')

    meanRate = meanRate (:,5);
    burstIndex = burstIndex (:,5);
    
    figure()
    scatter(meanRate,burstIndex)
        title('Burst Index Vs Mean Rate: Sleep Trials')
        xlabel('Mean Rate (seconds)')
        set(gca, 'xscale', 'log')
        ylabel('Burst Index')


% % 3D plot using the first trial for mean rate and burst index 
%     figure()
%     X = horzcat(meanRate(:, 1), burstIndex(:, 1), corecorded);
%     colours = hsv(5);
%     indices = kmeans(X, 4);
%     markerColours = zeros(length(X), 3);
%    
%     for jj=1: length(markerColours)
%         idx = indices(jj, 1);
%         if isnan(idx)
%             idx = 5;
%         end
%         markerColours(jj,1)=colours(idx, 1);
%         markerColours(jj,2)=colours(idx, 2);
%         markerColours(jj,3)=colours(idx, 3);
%     end
%     scatter3(meanRate(:, 1),burstIndex(:, 1), corecorded, 50, markerColours)
%         title('Burst Index Vs Mean Rate Vs corecorded cell count: All trials')
%         xlabel('Mean Rate (seconds)')
%         set(gca, 'xscale', 'log')
%         ylabel('Burst Index')
%         zlabel('Corecorded cell count')
% 
% % box plots from cluster data
%     figure()
%     boxplot(meanRate(:,2), indices)
%     ylabel('Mean Rate (seconds)')
%     set(gca, 'yscale', 'log')
%     
%     figure()
%     boxplot(burstIndex(:,2), indices)
%     ylabel('Burst Index')
%     
%     figure()
%     boxplot(corecorded, indices)
%     ylabel('Co-recorded cells')
%     
%     
% % box plots for mean rate, burst index and co-recorded cells
% 
%     figure()
%     boxplot(meanRate(:))
%     xlabel('All Trials')
%     ylabel('Mean Rate (seconds)')
%     set(gca, 'yscale', 'log')
% 
%     figure ()
%     boxplot(burstIndex(:))
%     xlabel('All Trials')
%     ylabel('Burst Index')
%     
%     figure ()
%     boxplot(corecorded)
%     xlabel('All Trials')
%     ylabel('Co-recorded cells')

% box plots for mean rate, burst index, and co-recorded cells per trial 
    
    
%     for jj = 1:5
%         figure()
%         boxplot(meanRate(:,jj), indices)
%         xlabel({'Clusters'},{':'}, string(env(1,jj)))
%         ylabel('Mean Rate (seconds)')
%         set(gca, 'yscale', 'log')
%     end
%     
%     for jj = 1:5
%         figure()
%         boxplot(burstIndex(:,jj),indices)
%         xlabel('Trials')
%         ylabel('Burst Index')

%         figure ()
%         boxplot(corecorded)
%         xlabel('All Trials')
% %         ylabel('Co-recorded cells')
%     end
    

    
end