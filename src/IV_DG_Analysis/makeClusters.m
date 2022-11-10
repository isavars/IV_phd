function [cluster1,cluster2,cluster3,cluster4] = makeClusters(spatData)
%CLUSTER makes clusters based on meanRate and burstIndex that can then be
%called by the other functions -- plan is to get this to use histology
%labels also - for now as a pre threshold but later it should be informing
%the values used for mean rate and burst index thresholds. The clusters
%should include CA3 cells also for now im only using data sets that were
%definitley not CA3.

%load ('allRats_spatData_CDB.mat', 'spatData')

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env; 
    SI_spat = spatData.SI_spat;
    nSpks = spatData.nSpks;
    
    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    nSpks(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;

    %find histology labels 

    cell_layer= histInfo();

    %make indexes for environment to use in comparisons with any file 
    maxFam = 2;
    FamInd = nan(length(meanRate),maxFam); 
    FamIndT = [];     
    NovInd = [];
    DiffInd = [];
    famCount = 0;
    numTrials = 5;

     for itCell= 1: length(meanRate)
        for itTrial = 1: numTrials
            if contains(cast(spatData.env(itCell,itTrial),'char'),'fam')
                FamIndT(itCell,itTrial) = itTrial;
                famCount = famCount + 1;
                if famCount <= maxFam
                     FamInd(itCell,famCount)= transpose(nonzeros(FamIndT(itCell,itTrial)));
                end
            elseif strcmp(cast(spatData.env(itCell,itTrial),'char'),'nov')
                NovInd(itCell,itTrial) = itTrial;
                NovInd = nonzeros(NovInd);
            elseif strcmp(cast(spatData.env(itCell,itTrial),'char'),'diff')
                DiffInd(itCell,itTrial) = itTrial;
                DiffInd = nonzeros(DiffInd);
            end 
        end
        famCount = 0;
     end

     %make rate change score and spatial information score
    
    for it_gma = 1: length (meanRate)
        awakeMeanRate(it_gma,:) = mean(meanRate(it_gma,(FamInd(1):FamInd(2))));
        SI_spat(it_gma,:) = max(SI_spat(it_gma,(1:end)));
    end

    sleepMeanRate = meanRate (:,end);
    
    rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));

    %make means of mean rate and burst index 

    meanMeanRate = [];
    meanBurstIndex = [];
    for it_gm = 1: length (meanRate)
        if nSpks(it_gm,FamInd(it_gm,1))  >= 10 || nSpks(it_gm,FamInd(it_gm,2))  >= 10 || nSpks(it_gm,NovInd(it_gm))  >= 10
            meanMeanRate(it_gm) = mean(meanRate(it_gm,end));% indexing here needs to generalize to all datasets 
            meanBurstIndex(it_gm) = mean(burstIndex(it_gm,end));
        end
    end


    %simple sorting based on histology label - this will hopefully produce
    %a distinct set of firing properties that can be used to sort the
    %dataset - from tetrodes and probes alike 

    granule_cluster = [];
    mossy_cluster = [];
    pyramidal_cluster = [];
    interneuron_cluster = []; %filter out super high firing rate from all the clusters to make this one 

    for itH = 1: length(cell_layer)
        if (meanMeanRate(itH) >= 1) && (meanBurstIndex(itH) <= 0.05) %interneurons - THIS LOOKS WRONG look up the knierim cuttoffs 
            interneuron_cluster = [interneuron_cluster;itH];
        else
            if strcmp(cell_layer{itH,1}, "GCL") 
               granule_cluster = [granule_cluster;itH];
            elseif strcmp(cell_layer{itH,1}, "HL")
                mossy_cluster = [mossy_cluster;itH];
            elseif strcmp(cell_layer{itH,1}, "CA3")
                pyramidal_cluster = [pyramidal_cluster;itH];
            else 
            end 
        end
    end

    cluster1 = interneuron_cluster;
    cluster2 = granule_cluster;
    cluster3 = mossy_cluster;
    cluster4 = pyramidal_cluster;

    
%     %make clusters based on firing properties 
% 
%     cluster1 = [];
%     cluster2 = [];
%     cluster3 = [];
%     cluster4 = [];
% 
%     for jj = 1: length(meanMeanRate)
%         if  (meanMeanRate(jj) >= 1.5) && (meanBurstIndex(jj) <= 0.05) %interneurons
%             cluster1 = [cluster1; jj];
%         elseif  (rateChange(jj) <= 0.3) && (meanMeanRate(jj) >= 0.01) && (meanMeanRate(jj) <= 1.5)  && (meanBurstIndex(jj) >= 0.01) && (meanBurstIndex(jj) < 0.15) %potential granule cells
%             cluster2 = [cluster2; jj];
%         elseif (meanMeanRate(jj) >= 0.02) && (meanBurstIndex(jj) >= 0.15) %potential mossy cells
%             cluster3 = [cluster3; jj];   
%         else 
%             cluster4 = [cluster4; jj];
%         end    
%     end 

end

