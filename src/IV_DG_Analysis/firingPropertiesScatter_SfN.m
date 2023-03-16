function firingPropertiesScatter_SfN (spatData) 

%gather data from spatData

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env; 
    SI_spat = spatData.SI_spat;
    SpkTs = spatData.SpkTs;
    
    meanRate(~(env=='fam' | env=='nov' | env=='diff' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov'| env=='diff'  | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='diff' | env=='sleep')) = NaN;
    
    %get burst index 
    for it_gm = 1: length (burstIndex)
        %meanMeanRate(it_gm) = mean(meanRate(it_gm,1:2));
        meanBurstIndex(it_gm) = mean(burstIndex(it_gm,1:2));
    end
%make indexes for environment to use in comparisons with any file 
    maxFam = 2;
    FamInd = nan(length(meanRate),maxFam); 
    FamIndT = [];     
    NovInd = [];
    famCount = 0;
    numTrials = 6;

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
            end 
        end
        famCount = 0;
     end
     FamInd1 = FamInd(:,1); %add these next time to simply a bit?
     FamInd2 = FamInd(:,2);  
     
%get clusters 


    [cluster1,cluster2,cluster3,cluster4] = makeClusters(spatData);
    

%make ? rate score
    
    for it_gma = 1: length (meanRate)
        meanMeanRate(it_gma,:) = mean(meanRate(it_gma,1:5));
    end

    awakeMeanRate = meanRate(:,1); %meanMeanRate;
    sleepMeanRate = meanRate (:,6);
    
    rateChange = sleepMeanRate ./ awakeMeanRate; %work on this tom said divide here ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));%
    
% make spatiality score

    for it_SI = 1: length (SI_spat)
        meanSI(it_SI,:) = mean(SI_spat(it_SI,1:3));
    end
    SI_spat = meanSI;

 %make special burst index from Senzai and Buzsaki

    burstIndex_Senzai = [];
    for itSP = 1: length(SpkTs)
        spike_AC = spk_crosscorr(cell2mat(SpkTs(itSP,FamInd(itSP,1))),'AC',0.001,0.3,900, 'norm', 'none');
        spike = mean(spike_AC(304:306));
        baseline= mean(spike_AC(501:601));
        burstIndex_Senzai = [burstIndex_Senzai; spike/baseline];
    end 

%make 2D scatter plots
    
%     figure()
%     scatter(meanMeanRate(cluster1), rateChange(cluster1), 'r')
%     hold on;
%         title('Sleep State Rate Change Vs Mean Rate')
%         xlabel('Mean Rate (seconds)')
%         set(gca, 'xscale', 'log')
%         ylabel('Sleep State Rate Change')
%     scatter(meanMeanRate(cluster2), rateChange(cluster2),'b')
%     hold on;
%     scatter(meanMeanRate(cluster3), rateChange(cluster3), 'g')
%     hold on;
%     scatter(meanMeanRate(cluster4), rateChange(cluster4), 'k')
%     hold off;
% 
% 
%     figure()
%     scatter(meanMeanRate(cluster1), meanBurstIndex(cluster1), 'r')
%     hold on;
%         title('Burst Index Vs Mean Rate')
%         xlabel('Mean Rate (seconds)')
%         set(gca, 'xscale', 'log')
%         ylabel('Burst Index')
%     scatter(meanMeanRate(cluster2), meanBurstIndex(cluster2),'b')
%     hold on;
%     scatter(meanMeanRate(cluster3), meanBurstIndex(cluster3), 'g')
%     hold on;
%     scatter(meanMeanRate(cluster4), meanBurstIndex(cluster4), 'k')
%     hold off;
% 
%     figure()
%     scatter(meanMeanRate(cluster1), SI_spat(cluster1), 'r')
%     hold on;
%         title('Spatiality Vs Mean Rate')
%         xlabel('Mean Rate (seconds)')
%         set(gca, 'xscale', 'log')
%         ylabel('Spatial Score')
%     scatter(meanMeanRate(cluster2), SI_spat(cluster2),'b')
%     hold on;
%     scatter(meanMeanRate(cluster3), SI_spat(cluster3), 'g')
%     hold on;
%     scatter(meanMeanRate(cluster4), SI_spat(cluster4), 'k')
%     hold off;
%     
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
    

%     figure()
%     scatter(burstIndex_Senzai(cluster1), rateChange(cluster1), 'r')
%     hold on;
%         title('Sleep State Rate Change Vs Burst Index')
%         xlabel('Burst Index')
%         ylabel('Sleep State Rate Change')
%         %set(gca, 'yscale', 'log')
%     scatter(burstIndex_Senzai(cluster2), rateChange(cluster2),'b')
%     hold on;
%     scatter(burstIndex_Senzai(cluster3), rateChange(cluster3), 'g')
%     hold on;
%     scatter(burstIndex_Senzai(cluster4), rateChange(cluster4), 'k')
%     hold off;
%     
    figure()
    scatter(meanBurstIndex(cluster1), rateChange(cluster1), 'r')
    hold on;
        title('Sleep State Rate Change Vs Burst Index')
        xlabel('Burst Index')
        ylabel('Sleep State Rate Change')
        set(gca, 'yscale', 'log')
    scatter(meanBurstIndex(cluster2), rateChange(cluster2),'b')
    hold on;
    scatter(meanBurstIndex(cluster3), rateChange(cluster3), 'g')
    hold on;
    scatter(meanBurstIndex(cluster4), rateChange(cluster4), 'k')
    hold off;
% %     
%     figure()
%     scatter(SI_spat(cluster1), rateChange(cluster1), 'r')
%     hold on;
%         title('Sleep State Rate Change Vs Spatiality')
%         xlabel('Spatial Score')
%         ylabel('Sleep State Rate Change')
%     scatter(SI_spat(cluster2), rateChange(cluster2),'b')
%     hold on;
%     scatter(SI_spat(cluster3), rateChange(cluster3), 'g')
%     hold on;
%     scatter(SI_spat(cluster4), rateChange(cluster4), 'k')
%     hold off;
%     
% % make 3D scatter plots 
%     
%     figure()
%     scatter3(meanMeanRate(cluster1),meanBurstIndex(cluster1),SI_spat(cluster1), 'r')
%         title('Burst Index Vs Mean Rate Vs Spatiality: All trials')
%         xlabel('Mean Rate (seconds)')
%         set(gca, 'xscale', 'log')
%         ylabel('Burst Index')
%         zlabel('Spatial Score')
%     hold on; 
%     scatter3(meanMeanRate(cluster2),meanBurstIndex(cluster2),SI_spat(cluster2), 'b')
%     hold on;
%     scatter3(meanMeanRate(cluster3),meanBurstIndex(cluster3),SI_spat(cluster3), 'g')
%     hold on;
%     scatter3(meanMeanRate(cluster4),meanBurstIndex(cluster4),SI_spat(cluster4), 'k')
%     hold off;
end 