function firingPropertiesScatter () 

load ('r889_All_spatData_fin.mat', 'spatData')

%gather data from spatData

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env; 
    SI_spat = spatData.SI_spat;
    
    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    
    for it_gm = 1: length (meanRate)
        meanMeanRate(it_gm) = mean(meanRate(it_gm,[1:3 5]));
        meanBurstIndex(it_gm) = mean(burstIndex(it_gm,[1:3 5]));
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
    cluster4 = nonzeros(cluster4);
    

%make ? rate score
    
    for it_gma = 1: length (meanRate)
        meanMeanRate(it_gma,:) = mean(meanRate(it_gma,1:3));
    end

    awakeMeanRate = meanMeanRate;
    sleepMeanRate = meanRate (:,5);
    
    rateChange = sleepMeanRate - awakeMeanRate; %work on this
    
% make spatiality score

    for it_SI = 1: length (SI_spat)
        meanSI(it_SI,:) = mean(SI_spat(it_SI,1:3));
    end
    SI_spat = meanSI;

%make 2D scatter plots
    
    figure()
    scatter(meanMeanRate(cluster1), rateChange(cluster1), 'r')
    hold on;
        title('Sleep State Rate Change Vs Mean Rate')
        xlabel('Mean Rate (seconds)')
        set(gca, 'xscale', 'log')
        ylabel('Sleep State Rate Change')
    scatter(meanMeanRate(cluster2), rateChange(cluster2),'b')
    hold on;
    scatter(meanMeanRate(cluster3), rateChange(cluster3), 'g')
    hold on;
    scatter(meanMeanRate(cluster4), rateChange(cluster4), 'k')
    hold off;

    figure()
    scatter(meanMeanRate(cluster1), SI_spat(cluster1), 'r')
    hold on;
        title('Spatiality Vs Mean Rate')
        xlabel('Mean Rate (seconds)')
        set(gca, 'xscale', 'log')
        ylabel('Spatial Score')
    scatter(meanMeanRate(cluster2), SI_spat(cluster2),'b')
    hold on;
    scatter(meanMeanRate(cluster3), SI_spat(cluster3), 'g')
    hold on;
    scatter(meanMeanRate(cluster4), SI_spat(cluster4), 'k')
    hold off;
    
    figure()
    scatter(meanBurstIndex(cluster1), SI_spat(cluster1), 'r')
    hold on;
        title('Spatiality Vs Burst Index')
        xlabel('Burst Index')
        ylabel('Spatial Score')
    scatter(meanBurstIndex(cluster2), SI_spat(cluster2),'b')
    hold on;
    scatter(meanBurstIndex(cluster3), SI_spat(cluster3), 'g')
    hold on;
    scatter(meanBurstIndex(cluster4), SI_spat(cluster4), 'k')
    hold off;
    
    figure()
    scatter(meanBurstIndex(cluster1), rateChange(cluster1), 'r')
    hold on;
        title('Sleep State Rate Change Vs Burst Index')
        xlabel('Burst Index')
        ylabel('Sleep State Rate Change')
    scatter(meanBurstIndex(cluster2), rateChange(cluster2),'b')
    hold on;
    scatter(meanBurstIndex(cluster3), rateChange(cluster3), 'g')
    hold on;
    scatter(meanBurstIndex(cluster4), rateChange(cluster4), 'k')
    hold off;
    
    figure()
    scatter(SI_spat(cluster1), rateChange(cluster1), 'r')
    hold on;
        title('Sleep State Rate Change Vs Spatiality')
        xlabel('Spatial Score')
        ylabel('Sleep State Rate Change')
    scatter(SI_spat(cluster2), rateChange(cluster2),'b')
    hold on;
    scatter(SI_spat(cluster3), rateChange(cluster3), 'g')
    hold on;
    scatter(SI_spat(cluster4), rateChange(cluster4), 'k')
    hold off;
    
% make 3D scatter plots 
    
    figure()
    scatter3(meanMeanRate(cluster1),meanBurstIndex(cluster1),SI_spat(cluster1), 'r')
        title('Burst Index Vs Mean Rate Vs Spatiality: All trials')
        xlabel('Mean Rate (seconds)')
        set(gca, 'xscale', 'log')
        ylabel('Burst Index')
        zlabel('Spatial Score')
    hold on; 
    scatter3(meanMeanRate(cluster2),meanBurstIndex(cluster2),SI_spat(cluster2), 'b')
    hold on;
    scatter3(meanMeanRate(cluster3),meanBurstIndex(cluster3),SI_spat(cluster3), 'g')
    hold on;
    scatter3(meanMeanRate(cluster4),meanBurstIndex(cluster4),SI_spat(cluster4), 'k')
    hold off;
end 