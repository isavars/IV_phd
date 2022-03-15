function [cluster1,cluster2,cluster3,cluster4] = makeClusters(meanRate, burstIndex)
%CLUSTER makes clusters based on meanRate and burstIndex

load ('allRatsDG_P18P23_temp.mat', 'spatData')

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env; 
%     SI_spat = spatData.SI_spat;
%     nSpks = spatData.nSpks;
    
    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
%     SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
%     nSpks(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;

    if nargin == meanRate && nargin == burstIndex
        
        cluster1 = [];
        cluster2 = [];
        cluster3 = [];
        cluster4 = [];

        for jj = 1: length(meanMeanRate)
            if  (meanMeanRate(jj) > 1) && (meanMeanRate(jj) < 10) && (meanBurstIndex(jj) < 0.05) %mystery cells
                cluster1 = [cluster1; jj];
            elseif  (meanMeanRate(jj) > 0.01) && (meanMeanRate(jj) < 0.99)  && (meanBurstIndex(jj) > 0.02) && (meanBurstIndex(jj) < 0.15) %potential granule cells
                cluster2 = [cluster2; jj];
            elseif (meanMeanRate(jj) > 0.02) && (meanBurstIndex(jj) > 0.15) %potential mossy cells
                cluster3 = [cluster3; jj];
            else 
                cluster4 = [cluster4; jj];
            end    
        end 
    end
end

