function waveformPCA_TW (spatData)

% 
%gather data from spatData

    meanRate = spatData.meanRate;
    env = spatData.env; 
    wf_means = spatData.wf_means;
    waveforms = spatData.waveforms;
    nSpks = spatData.nSpks;
    SpkTs = spatData.SpkTs;
    
    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    
%idxs for environment to use in comparisons with any file 
    maxFam = 2;
    FamInd = nan(length(meanRate),maxFam); 
    FamIndT = [];     
    NovInd = [];
    diffInd =[];
    famCount = 0;
    trialNum = 6;

    for itCell= 1: length(meanRate)
        for itTrial = 1: length(trialNum)
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
                diffInd(itCell,itTrial) = itTrial;
                diffInd = nonzeros(diffInd);
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

%make spike width  

    for itSp = 1: length (nSpks) 
        [~, maxSpksPos] = max(nSpks(itSp,1:3));
        WFs (itSp,:) = wf_means(itSp, maxSpksPos);
    end

    InCluster = [];
    ExCluster = [];
    for itWF = 1: length (WFs) 
        if length(cell2mat(waveforms(itWF,FamInd(itWF,1)))) == 50
            spikeFeatures = spk_characterisewaveform(WFs{itWF,:}, [100,100,100,100]);%save this somewhere else it takes ages 
            spikeWidth (itWF, 1) = spikeFeatures.pk2TrTime;
            spikeWidth (itWF,2) = spikeFeatures.pk2TrAmp;
        else 
            spikeWidth (itWF, 1) = NaN;
            spikeWidth (itWF,2) = NaN;
        end
        if spikeWidth (itWF,1) < 0.425 && burstIndex_Senzai(itWF) <= 1
            InCluster = [InCluster;itWF]; 
        else 
            ExCluster = [ExCluster;itWF];
        end
    end
    
%PCA on excitatory cells only

    waveform_inds = sub2ind(size(waveforms), ExCluster,FamInd(ExCluster,1));
    ex_waveforms = waveforms(waveform_inds);

    pca_data = [];
    for itEx = 1: length(ex_waveforms)
         if ~isnan(cell2mat(ex_waveforms(itEx))) & length(cell2mat(ex_waveforms(itEx))) == 50
             pca_data = [pca_data; interp1(1:50, ex_waveforms{itEx}(1:50),1:0.48:50,'spline')];
         end
    end
    
    % Normalise so each WF peak is 1.
    pca_data      = pca_data ./ max(pca_data,[],2);
    
    % Shift WFs so that peaks aligned.
    [~,pkInd]  = max(pca_data,[],2);
    firstPkInd = min(pkInd);
    nSampsWF   = size(pca_data,2);
    padWF      = size(pca_data);
    for itWF=1:size(pca_data,1)
        WFShift = pkInd(itWF) - firstPkInd;
        padWF(itWF, 1:(nSampsWF-WFShift) ) = pca_data( itWF, (WFShift+1):nSampsWF );
    end
    % Use whole WF, or only from peak onwards? (Some parts of S&B imply the latter). Seems worse the whole WF, though.
    if 1
        pca_data = padWF( : , 1:round(nSampsWF*0.8) );  % Cut WF directly to 0.8ms here, should remove NaNs.
    else
        pca_data = padWF( : , firstPkInd:round(nSampsWF*0.8) );  % Cut WF directly to 0.8ms here, should remove NaNs.
    end
    
    % Take 2nd differential
    diff_pca_data = diff(pca_data,2,2);
%     for itPCA = 1: size(pca_data,1)
%         diff_pca_data = [diff_pca_data; diff(pca_data(itPCA,:),2)];
%     end
    
    % Mean normalise so every WF mean is zero (doesn't seem to mak much difference)
    diff_pca_data = diff_pca_data - mean( diff_pca_data, 2 );
    

    coeff =  pca(diff_pca_data);
    w_PC1 = coeff(:,1);
    w_PC2 = coeff(:,2);
    % Plot PC weights.
    figure;
    for itPC = 1:2
        subplot(1,2,itPC);   
        plot( coeff(:,itPC) );   
        hold on;
        plot( [1 1].*firstPkInd, [-0.4 0.4], 'k:' );
    end

    %dot product of 2nd derrivative of upsampled waveforms with w-PC1 and w-PC2.
    % Get dot product 'by hand' - is just each element in vector A times each in vector B, summed together.
    wfPC1 = sum( diff_pca_data.*(w_PC1'), 2 );
    wfPC2 = sum( diff_pca_data.*(w_PC2'), 2 );
%     for itPD = 1: size(diff_pca_data,1)
%          wfPC1 = [wfPC1; dot(diff_pca_data(itPD,:).',w_PC1)];
%          wfPC2 = [wfPC2; dot(diff_pca_data(itPD,:).',w_PC2)];
%     end

    % Plot weight distributions.
    figure;
    subplot(2,2,1);  histogram(wfPC1)
    subplot(2,2,2);  histogram(wfPC2)
    sAx = subplot(2,2,3);  
    scatter(wfPC1,wfPC2);
    xEdge = linspace( sAx.XLim(1), sAx.YLim(2), 40 );
    yEdge = linspace( sAx.YLim(1), sAx.YLim(2), 40 );
    cnts  = histcounts2( wfPC1,wfPC2, xEdge, yEdge );
    subplot(2,2,4);  imagesc(cnts');  
    set(gca, 'ydir', 'normal');  %   'xticklabel', xEdge(get(gca,'xtick')), 'yticklabel', yEdge(get(gca,'ytick'))
    
    % Divide into groups and plot WFs.
    grInd{1} = wfPC1>-0.05; % & wfPC1<-0.03;
    grInd{2} = wfPC1<-0.05;

%     grInd{1} = wfPC2>-0.015;
%     grInd{2} = wfPC2<-0.015;


    figure;
    for itGr = 1:length(grInd)
        plot( mean(pca_data( grInd{itGr}, : ) ) );
        hold on;
    end
    
    

end