function class_waveformPCA (data)
% makes waveform PCA features from Senzai and Buzsaki 2017
% INPUT - spatData 

load (data, 'spatData')

%gather data from spatData

    meanRate = spatData.meanRate;
    env = spatData.env; 
    wf_means = spatData.wf_means;
    waveforms = spatData.waveforms;
    nSpks = spatData.nSpks;
    SpkTs = spatData.SpkTs;
    TP_latency = mean(spatData.TP_latency,2);
    
%     meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN; % not sure
%     why this line is here 
    
%idxs for environment to use in comparisons with any file - are these
%necessary? 

    maxFam = 2; %i guess some might be missing a third fam trial? 
    FamInd = nan(length(meanRate),maxFam); 
    FamIndT = [];     
    NovInd = [];
    famCount = 0;

    for itCell= 1: length(meanRate)
        for itTrial = 1: size(env,2) 
            if contains(cast(env(itCell,itTrial),'char'),'fam')
                FamIndT(itCell,itTrial) = itTrial;
                famCount = famCount + 1;
                if famCount <= maxFam
                     FamInd(itCell,famCount)= transpose(nonzeros(FamIndT(itCell,itTrial)));
                end
            elseif strcmp(cast(env(itCell,itTrial),'char'),'nov')
                NovInd(itCell,itTrial) = itTrial;
                NovInd = nonzeros(NovInd);
            end 
        end
        famCount = 0;
    end
    
%make special burst index from Senzai and Buzsaki

    burstIndex_Senzai = [];
    for itSP = 1: length(SpkTs)
        spike_AC = spk_crosscorr(cell2mat(SpkTs(itSP,FamInd(itSP,1))),'AC',0.001,0.3,900, 'norm', 'none');
        spike = mean(spike_AC(304:306));
        baseline= mean(spike_AC(501:601));
        burstIndex_Senzai = [burstIndex_Senzai; spike/baseline];
    end 

%make spike width - this is trough-to-peak latency which is now on spatData
%can change if needed. 

%chose waveform with most spikes 
    for itSp = 1: length (nSpks) 
        [~, maxSpksPos] = max(nSpks(itSp,1:3));
        WFs (itSp,:) = wf_means(itSp, maxSpksPos);
    end

%filter out interneurons 
    InCluster = [];
    ExCluster = [];
%toms method 
%     for itWF = 1: length (WFs) 
%         if length(cell2mat(waveforms(itWF,FamInd(itWF,1)))) == 50
%             spikeFeatures = spk_characterisewaveform(WFs{itWF,:}, [100,100,100,100]);%save this somewhere else it takes ages 
%             spikeWidth (itWF, 1) = spikeFeatures.pk2TrTime;
%             spikeWidth (itWF,2) = spikeFeatures.pk2TrAmp;
%         else 
%             spikeWidth (itWF, 1) = NaN;
%             spikeWidth (itWF,2) = NaN;
%         end
%         if spikeWidth (itWF,1) < 0.425 && burstIndex_Senzai(itWF) <= 1.2
%             InCluster = [InCluster;itWF]; 
%         else 
%             ExCluster = [ExCluster;itWF]; 
%         end
%     end
%using my trough to peak (max is peak and min is trough of the wf and its a
%subtraction in ms. 
    for itWF = 1: length (WFs)
        if TP_latency(itWF) < 0.425 && burstIndex_Senzai(itWF) <= 1.2
            InCluster = [InCluster;itWF]; 
        else 
            ExCluster = [ExCluster;itWF]; 
        end            
    end 
    
% PCA on excitatory cells only.Paper says - Waveform PCA was performed using the time between 
    % 0 to 0.8ms period of the second derivative of the up-sampled average 
    % waveform. 

    waveform_inds = sub2ind(size(waveforms), ExCluster,FamInd(ExCluster,1));
    ex_waveforms = waveforms(waveform_inds);
    pca_data = [];
    for itEx = 1: length(ex_waveforms)
         if ~isnan(cell2mat(ex_waveforms(itEx))) % && length(cell2mat(ex_waveforms(itEx))) == 97 %changed from 50 don't know why this is here 
             pca_data = [pca_data; interp1(1:97, cell2mat(ex_waveforms(itEx)),1:0.48:97,'spline')];%[pca_data; interp1(1:50, cell2mat(ex_waveforms(itEx)),1:0.48:50,'spline')]
         end
    end
    diff_pca_data = [];
    for itPCA = 1: size(pca_data,1)
        diff_pca_data = [diff_pca_data; diff(pca_data(itPCA,20:100),2)]; %inserted 1:80 for the .8 ms time window
    end
    coeff =  pca(diff_pca_data);
    w_PC1 = coeff(:,1);
    w_PC2 = coeff(:,2);

    %dot product of 2nd derrivative of upsampled waveforms with w-PC1 and
    %w-PC2.
 
    wfPC1 = [];
    wfPC2 = [];
    for itPD = 1: size(diff_pca_data,1)
         wfPC1 = [wfPC1; dot(diff_pca_data(itPD,:).',w_PC1)];
         wfPC2 = [wfPC2; dot(diff_pca_data(itPD,:).',w_PC2)];
    end
    figure()
    plot(w_PC1)
    xlabel('Samples')
    ylabel('w-PC1')
    title('Principal Component 1')
    
    figure()
    histogram(wfPC1)
    xlabel('w-PC1')
    ylabel('Cell Count')
    title('Histogram of PC1')
    
    figure()
    plot(w_PC2)
    xlabel('Samples')
    ylabel('w-PC2')
    title('Principal Component 2')
    
    figure()
    histogram(wfPC2)
    xlabel('w-PC2')
    ylabel('Cell Count')
    title('Histogram of PC2')
    
    figure()
    scatter(wfPC1,wfPC2)
    xlabel('w-PC1')
    ylabel('w-PC2')
    title('Scatter plot of PC1 vs PC2')

    

end