function firingPropertiesScatter (data) 

load (data, 'spatData')

%gather data from spatData

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env; 
    SI_spat = spatData.SI_spat;
    wf_means = spatData.wf_means;
    waveforms = spatData.waveforms;
    nSpks = spatData.nSpks;
    SpkTs = spatData.SpkTs;
    
    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;



%make indexes for environment to use in comparisons with any file 
    maxFam = 3;%iv changed to 3 from 2 02/12/21
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
%make rate change score
    
    for it_gma = 1: length (meanRate)
        awakeMeanRate(it_gma,:) = mean(meanRate(it_gma,(FamInd(it_gma,1):FamInd(it_gma,2))));
    end

    sleepMeanRate = meanRate (:,5);
    
    rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));

%gather age data from cellInfo 

    cellInfo = getCellInfo();
    corecorded = cellInfo(:,4);
    
%make clusters 

    meanMeanRate = [];
    meanBurstIndex = [];
    for it_gm = 1: length (meanRate)
        if nSpks(it_gm,FamInd(it_gm,1))  >= 10 || nSpks(it_gm,FamInd(it_gm,2))  >= 10 || nSpks(it_gm,NovInd(it_gm))  >= 10 % from final rateOverlap file
            meanMeanRate(it_gm) = mean(meanRate(it_gm,5));% indexing here needs to generalize to all datasets 
            meanBurstIndex(it_gm) = mean(burstIndex(it_gm,5));
        end
    end

    cluster1 = [];
    cluster2 = [];
    cluster3 = [];
    cluster4 = [];
    cluster5 = [];

    for jj = 1: length(meanMeanRate)
%         if  (meanMeanRate(jj) >= 10)   %interneurons (meanMeanRate(jj) >= 1) && (meanMeanRate(jj) <= 10) && (meanBurstIndex(jj) <= 0.05)
%             cluster1 = [cluster1; jj];
%         elseif  (rateChange(jj) <= 0.3) && (meanMeanRate(jj) >= 0.01) && (meanMeanRate(jj) <= 1.5)  && (meanBurstIndex(jj) >= 0.01) && (meanBurstIndex(jj) < 0.15) %potential granule cells burst index min 0.02 (from rate overlap function)
%             cluster2 = [cluster2; jj];
%         elseif (meanMeanRate(jj) <= 1.5) && (meanBurstIndex(jj) >= 0.15) && (corecorded(jj) > 3) %potential CA3 cells
%             cluster4 = [cluster4; jj];
%         elseif (meanMeanRate(jj) >= 0.02) && (meanBurstIndex(jj) >= 0.15) %potential mossy cells
%             cluster3 = [cluster3; jj];
%         else 
%             cluster5 = [cluster5; jj];
%         end    
            if  (meanMeanRate(jj) > 1)  && (meanBurstIndex(jj) < 0.05) || (meanMeanRate(jj) > 10)%interneurons && (meanMeanRate(jj) < 10)
                cluster1 = [cluster1; jj];
            elseif  (rateChange(jj) <= 0.3) && (meanMeanRate(jj) > 0.01) && (meanMeanRate(jj) < 1.5)  && (meanBurstIndex(jj) > 0.02) && (meanBurstIndex(jj) < 0.15) %potential granule cells
                cluster2 = [cluster2; jj];
            elseif (meanMeanRate(jj) > 0.02) && (meanBurstIndex(jj) > 0.15) && (meanMeanRate(jj) < 10)%potential mossy cells
                cluster3 = [cluster3; jj];
            else 
                cluster4 = [cluster4; jj];
            end
    end 
     

%make rate change score
    
    for it_gma = 1: length (meanRate)
        awakeMeanRate(it_gma,:) = mean(meanRate(it_gma,(FamInd(it_gma,1):FamInd(it_gma,2))));
    end

    sleepMeanRate = meanRate (:,5);
    
    rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));%awakeMeanRate-sleepMeanRate; %
   
    % removing outliers 
%     for itRC = 1: length(rateChange)
%         if rateChange(itRC) <= -1000
%             rateChange(itRC) = [];
%         end
%     end
%     
% make spatiality score

    for it_SI = 1: length (meanRate)
        meanSI(it_SI,:) = mean(SI_spat(it_SI,2:4));
    end
    SI_spat = meanSI;
    
%make special burst index from Senzai and Buzsaki

%     burstIndex_Senzai = [];
%     for itSP = 1: length(SpkTs)
%         spike_AC = spk_crosscorr(cell2mat(SpkTs(itSP,FamInd(itSP,1))),'AC',0.001,0.3,900, 'norm', 'none');
%         spike = mean(spike_AC(304:306));
%         baseline= mean(spike_AC(501:601));
%         burstIndex_Senzai = [burstIndex_Senzai; spike/baseline];
%     end 

%make spike width  

    for itSp = 1: length (nSpks) 
        [~, maxSpksPos] = max(nSpks(itSp,1:3));
        WFs (itSp,:) = wf_means(itSp, maxSpksPos);
    end

%     InCluster = [];
%     ExCluster = [];
%     for itWF = 1: length (WFs) 
%         if length(cell2mat(waveforms(itWF,FamInd(itWF,1)))) == 50
%             spikeFeatures = spk_characterisewaveform(WFs{itWF,:}, [100,100,100,100]);%save this somewhere it takes ages 
%             spikeWidth (itWF, 1) = spikeFeatures.pk2TrTime;
%             spikeWidth (itWF,2) = spikeFeatures.pk2TrAmp;
%         else 
%             spikeWidth (itWF, 1) = NaN;
%             spikeWidth (itWF,2) = NaN;
%         end
%         if spikeWidth (itWF,1) < 0.425 && burstIndex_Senzai(itWF) <= 1
%             InCluster = [InCluster;itWF]; 
%         else 
%             ExCluster = [ExCluster;itWF]; 
%         end
%     end
    
% PCA on excitatory cells only

%     waveform_inds = sub2ind(size(waveforms), ExCluster,FamInd(ExCluster,1));
%     ex_waveforms = waveforms(waveform_inds);
%     pca_data = [];
%     for itEx = 1: length(ex_waveforms)
%          if ~isnan(cell2mat(ex_waveforms(itEx))) & length(cell2mat(ex_waveforms(itEx))) == 50
%              pca_data = [pca_data; interp1(1:50, cell2mat(ex_waveforms(itEx)),1:0.48:50,'spline')];
%          end
%     end
%     diff_pca_data = [];
%     for itPCA = 1: size(pca_data,1)
%         diff_pca_data = [diff_pca_data; diff(pca_data(itPCA,:),2)];
%     end
%     figure()
%     plot(diff_pca_data(1,:))
%     coeff =  pca(diff_pca_data);
%     w_PC1 = coeff(:,1);
%     figure()
%     plot(w_PC1)
%     
%     w_PC2 = coeff(:,2);
%     figure()
%     plot(w_PC2)
    %dot product of 2nd derrivative of upsampled waveforms with w-PC1 and
    %w-PC2    
 
%     wfPC1 = [];
%     wfPC2 = [];
%     for itPD = 1: size(diff_pca_data,1)
%          wfPC1 = [wfPC1; dot(diff_pca_data(itPD,:).',w_PC1)];
%          wfPC2 = [wfPC2; dot(diff_pca_data(itPD,:).',w_PC2)];
%     end
%     figure()
%     histogram(wfPC1)
%     figure()
%     histogram(wfPC2)
%     figure()
%     scatter(wfPC1,wfPC2)
    
    
    
%spike width vs mean rate scatter 

%     figure()
%     scatter(meanRate(:,FamInd(1)), spikeWidth(:,2))
%         title('Spike Width Vs Mean Rate')
%         xlabel('Mean Firing Rate (Hz)')
%         set(gca, 'xscale', 'log')
%         ylabel('Spike width (ms)')
    
%     figure()
%     scatter(spikeWidth(:,1),burstIndex_Senzai)%scatter(spikeWidth(:,1), burstIndex(:,FamInd(1)))%
%         title('Spike width Vs Burst Index')
%         ylabel('Burst Index')
%         set(gca, 'yscale', 'log')
%         xlabel('Spike width (ms)')
%         
%     figure()
%     histogram(burstIndex_Senzai,10)
        
%     figure()
%     scatter(burstIndex(:,FamInd(1)), spikeWidth(:,2))
%         title('Peak to Trough Amplitide Vs Burst Index')
%         xlabel('Burst Index')
%         ylabel('Peak to Trough Amplitude')
    
    
% rate change histogram 
    
%     figure()
%     hold all;
%     h1 = histfit(rateChange(cluster1),10, 'kernel');
% %         axis ([-200 200 0 5])
%         title('State Dependent Rate Change Distribution')
%         xlabel('State Dependent Rate Change')
%         ylabel('Frequency')
%     h1(2).Color = 'r';
%     delete (h1(1)); 
%     h2 = histfit(rateChange(cluster2),10, 'kernel');
%     delete (h2(1));
%     h2(2).Color = 'b';
%     h3 = histfit(rateChange(cluster3),10, 'kernel');
%     delete (h3(1));
%     h3(2).Color = 'g';
%     h4 = histfit(rateChange(cluster4),10, 'kernel');
%     delete (h4(1));
%     h4(2).Color = 'k';
%     legend('Putative Interneurons','Putative Granule','Putative Mossy', 'Unclassified');
    
%make 2D scatter plots

%     figure()
%     scatter(meanMeanRate(cluster1), rateChange(cluster1), 'r')
%     hold all;
%         title('State Dependent Rate Change Vs Mean Rate')
%         xlabel('Mean Rate (seconds)')
%         set(gca, 'xscale', 'log')
%         ylabel('State Dependent Rate Change')
%     scatter(meanMeanRate(cluster2), rateChange(cluster2),'b')
%     scatter(meanMeanRate(cluster3), rateChange(cluster3), 'g')
%     scatter(meanMeanRate(cluster4), rateChange(cluster4), 'k')
%     legend('Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');

    figure()
    
    hold on;
        title('Burst Index Vs Mean Rate')
        xlabel('Mean Rate (seconds)')
        set(gca, 'xscale', 'log')
        ylabel('Burst Index')
%     scatter(meanMeanRate(InCluster), meanBurstIndex(InCluster), 'r')
%     scatter(meanMeanRate(ExCluster), meanBurstIndex(ExCluster), 'b')
    scatter(meanMeanRate(cluster1), meanBurstIndex(cluster1),'r', 'LineWidth',1.5)
    scatter(meanMeanRate(cluster2), meanBurstIndex(cluster2), 'MarkerEdgeColor',[0.9290 0.6940 0.1250],'LineWidth',1.5)
    scatter(meanMeanRate(cluster3), meanBurstIndex(cluster3), 'MarkerEdgeColor',[0.4660 0.6740 0.1880],'LineWidth',1.5)
    scatter(meanMeanRate(cluster4), meanBurstIndex(cluster4), 'k')
    scatter(meanMeanRate(cluster5), meanBurstIndex(cluster5),'k')
%     legend('Putative Inhibitory','Putative Exitatory');
    legend('Putative Interneurons','Putative Granule','Putative Mossy', 'Unclassified');
    
    
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
%     legend('Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
%     hold off;
%     
%     figure()
%     scatter(meanBurstIndex(cluster1), SI_spat(cluster1), 'r')
%     hold on;
%         title('Spatiality Vs Burst Index')
%         xlabel('Burst Index')
%         ylabel('Spatial Score')
%     scatter(meanBurstIndex(cluster2), SI_spat(cluster2),'b')
%     hold on;
%     scatter(meanBurstIndex(cluster3), SI_spat(cluster3), 'g')
%     hold on;
%     scatter(meanBurstIndex(cluster4), SI_spat(cluster4), 'k')
%     legend('Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
%     hold off;
    
%     figure()
%     scatter(meanBurstIndex(cluster1), rateChange(cluster1), 'r')
%     hold on;
%         title('State Dependent Rate Change Vs Burst Index')
%         xlabel('Burst Index')
%         ylabel('State Dependent Rate Change')
%     scatter(meanBurstIndex(cluster2), rateChange(cluster2),'b')
%     hold on;
%     scatter(meanBurstIndex(cluster3), rateChange(cluster3), 'g')
%     hold on;
%     scatter(meanBurstIndex(cluster4), rateChange(cluster4), 'k')
%     legend('Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
%     hold off;
%     
%     figure()
%     scatter(SI_spat(cluster1), rateChange(cluster1), 'r')
%     hold on;
%         title('State Dependent Rate Change Vs Spatiality')
%         xlabel('Spatial Score')
%         ylabel('State Dependent Rate Change')
%         
%     scatter(SI_spat(cluster2), rateChange(cluster2),'b')
%     hold on;
%     scatter(SI_spat(cluster3), rateChange(cluster3), 'g')
%     hold on;
%     scatter(SI_spat(cluster4), rateChange(cluster4), 'k')
%     legend('Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
%     hold off;
    
% make 3D scatter plots 
    
%     figure()
%     scatter3(meanMeanRate(cluster1),meanBurstIndex(cluster1),SI_spat(cluster1), 'r')
%         title('Burst Index Vs Mean Rate Vs Spatiality')
%         xlabel('Mean Rate (seconds)')
%         set(gca, 'xscale', 'log')
%         ylabel('Burst Index')
%         zlabel('Spatial Score')
%     hold all; 
%     scatter3(meanMeanRate(cluster2),meanBurstIndex(cluster2),SI_spat(cluster2), 'b')
%     scatter3(meanMeanRate(cluster3),meanBurstIndex(cluster3),SI_spat(cluster3), 'g')
%     scatter3(meanMeanRate(cluster4),meanBurstIndex(cluster4),SI_spat(cluster4), 'k')
%     legend('Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
%     
%     figure()
%     scatter3(meanMeanRate(cluster1),meanBurstIndex(cluster1),rateChange(cluster1), 'r')
%         title('Burst Index Vs Mean Rate Vs State Dependent Rate Change')
%         xlabel('Mean Rate (seconds)')
%         set(gca, 'xscale', 'log')
%         ylabel('Burst Index')
%         zlabel('State Dependent Rate Change')
%     hold all; 
%     scatter3(meanMeanRate(cluster2),meanBurstIndex(cluster2),rateChange(cluster2), 'b')
%     scatter3(meanMeanRate(cluster3),meanBurstIndex(cluster3),rateChange(cluster3), 'g')
%     scatter3(meanMeanRate(cluster4),meanBurstIndex(cluster4),rateChange(cluster4), 'k')
%     legend('Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
%     
%     figure()
%     scatter3(meanMeanRate(cluster1),SI_spat(cluster1),rateChange(cluster1), 'r')
%         title('Spatiality Vs Mean Rate Vs State Dependent Rate Change')
%         xlabel('Mean Rate (seconds)')
%         set(gca, 'xscale', 'log')
%         ylabel('Spatial Score')
%         zlabel('State Dependent Rate Change')
%     hold all; 
%     scatter3(meanMeanRate(cluster2),SI_spat(cluster2),rateChange(cluster2), 'b')
%     scatter3(meanMeanRate(cluster3),SI_spat(cluster3),rateChange(cluster3), 'g')
%     scatter3(meanMeanRate(cluster4),SI_spat(cluster4),rateChange(cluster4), 'k')
%     legend('Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
%     
%     figure()
%     scatter3(SI_spat(cluster1),meanBurstIndex(cluster1),rateChange(cluster1), 'r')
%         title('Burst Index Vs Spatiality Vs State Dependent Rate Change')
%         xlabel('Spatial Score')
%         ylabel('Burst Index')
%         zlabel('State Dependent Rate Change')
%     hold all; 
%     scatter3(SI_spat(cluster2),meanBurstIndex(cluster2),rateChange(cluster2), 'b')
%     scatter3(SI_spat(cluster3),meanBurstIndex(cluster3),rateChange(cluster3), 'g')
%     scatter3(SI_spat(cluster4),meanBurstIndex(cluster4),rateChange(cluster4), 'k')
%     legend('Mystery Cells','Putative Granule','Putative Mossy', 'Unclassified');
  
    
end 