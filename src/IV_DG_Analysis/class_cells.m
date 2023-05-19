function class_cells(data,electrodes)
    %this function needs to read in cells from spatData and electrode
    %position table (which has DS2 info and other postion estimate methods)
    %and provide cell classification ids (the length of spatData) for 
    %interneuron, granule cell, mossy cell and CA3 pyramidal cell. 

    % ISSUE #1 - might need to load tetrode index from dat2spikes to have a
    % way to interpret what channels belong to what tetrode 
    % TO DO - Add linear probe infon and distance from inversion info to
    %         DS2_info.
    %       - sort out how the order works where the clusters still reference
    %         spatData

    %load spatial Data 
    load (data, 'spatData');
    %load useful parts from spatData
    meanRate = spatData.meanRate;
    animal = spatData.animal;
    wf_means = spatData.wf_means;
    waveforms = spatData.waveforms; %waveforms is created from the mean wf with the maximum amplitude - need to keep channel id for it to use DS2 labels
    max_wf_chan = spatData.max_wf_channel; %this is per trial needs to be one value 
    nSpks = spatData.nSpks;
    SpkTs = spatData.SpkTs;
    TP_latency = mean(spatData.TP_latency,2);
    %load electrode Positions 
    load (electrodes, 'elePos');
    rat_id = cellstr(elePos.rat_id);
    hist_labels = elePos.hist_labels;

    %chose trial with most spikes to use wfs, spike times, mean firing rate, max wf channel from it 
    for itSp = 1: length (nSpks) 
        [~, maxSpksPos] = max(nSpks(itSp,1:5)); % all wake trials (hc change somewhere)
        WFs (itSp,:) = wf_means(itSp, maxSpksPos); %gets best wf from wf means
        max_wf_chan(itSp,:) = max_wf_chan(itSp,maxSpksPos); 
        burstIndex (itSp,:) = spatData.burstIndex(itSp, maxSpksPos);
        STs (itSp,:) = SpkTs(itSp, maxSpksPos); %gets the spiketimes from the same trial to calculate the burst index
%         awakeMeanRate (itSp,:) = meanRate(itSp, maxSpksPos); %gets the wake mean firing rate for the most active wake trial (doesn this make sense?)
    end
    
    %get tetrode to shank and channel mapping to be used by cluster feature
    %functions 
    [tetShankChan] = makeTetShankChan(spatData,max_wf_chan);


    %step 1 - filter by gross histology postion CA3 vs DG 
        % reads elePos values for hist_labels per tetrode and if a cell is
        % in the DG it gets labeled DG_cluster and if its in CA3
        % CA3_cluster. 
    
     DG_cluster = [];
     CA3_cluster = []; 

     for it_cells = 1: height(spatData)
         for it_ep = 1: height(elePos)
             if strcmp(animal(it_cells),rat_id(it_ep))
                 hist_label = hist_labels(it_ep,tetShankChan(it_cells,2));
                 if strcmp(hist_label, "DG")
                    DG_cluster = [DG_cluster;it_cells];
                 elseif strcmp(hist_label, "CA3")
                    CA3_cluster = [CA3_cluster;it_cells]; 
                 end
             end
         end
     end  
    % this is a branching point for DG and CA3 - the CA3 then needs to be
    % sorted into interneurons and CA3 pyramidal cells also - should It all
    % be in one big loop that iterates over rows of spatData? 


    %step 2 - filter by excitatory vs inhibitory using trough-to-peak measure

        %make burst index from Senzai and Buzsaki 
        burstIndex_Senzai = [];
        for itSP = 1: length(SpkTs)
            spike_AC = spk_crosscorr(cell2mat(STs(itSP)),'AC',0.001,0.3,900, 'norm', 'none'); 
            spike = mean(spike_AC(304:306));
            baseline= mean(spike_AC(501:601));
            burstIndex_Senzai = [burstIndex_Senzai; spike/baseline];
        end     
        %make the inhibitory and excitatory clusters for CA3 and DG subgroups
        %-mantaining row numbers from spatData
        InCluster = [];
        ExCluster = [];
        for itWF = 1: length (WFs)
            if TP_latency(itWF) < 0.425 && burstIndex_Senzai(itWF) <= 1.2
                InCluster = [InCluster;itWF]; 
            else 
                ExCluster = [ExCluster;itWF]; 
            end            
        end 
        % temporarily making a DG_Ex_cluster
        DG_ExCluster = DG_cluster(ismember(DG_cluster, ExCluster));

    %step 3 - call in position relative to DS2 (specifically DS2 amplitude
    %for the ones that have it and more informed histology lables for the
    %ones that dont) 
        %loop through spatData(DG_ExCluster) and find the DS2 amplitudes of the chosen
        %cells (can select per channel using wf_max_chan index) 

        DS2_amplitudes = [];
        for it_DG_Ex = DG_ExCluster'
            for it_ep = 1: height(elePos)
                if strcmp(animal(it_DG_Ex),rat_id(it_ep))
                    DS2_amplitude = elePos.DS2_amplitude(it_ep,tetShankChan(it_DG_Ex,3)); %third column of tetShankChan is for channel
                    DS2_amplitudes = [DS2_amplitudes;DS2_amplitude];
                end
            end 
        end

    %step 4 - create wf_PCA features (call waveformPCA, run as a
    %subfunciton) need to know what outputs from the pca are the right ones
    %to use here. 

        [wfPC1,wfPC2] = waveformPCA(DG_ExCluster,waveforms);

    %step 5 - create sleep vs wake firing rate (rateChange) - needs to be
    %indexed within the cluster for DG or CA3 and withought Interneurons -
    %need to make a DG_exCluster and a CA3_exCluster - or do this
    %differently with one big loop 

        awakeMeanRate = nanmean(meanRate(:,1:5),2);
        awakeMeanRate = awakeMeanRate(DG_ExCluster);
        sleepMeanRate = meanRate (DG_ExCluster,end);
        rateChange = awakeMeanRate ./ sleepMeanRate;

    %step 6 - run a second PCA on the 4 features wfPC1, wfPC2, rate change
    %and DS2 amplitude

    % step 7 -run k means as suggested by S&B - once the correct PCA
    %features are identified. 
    
    burstIndex = burstIndex(DG_ExCluster);

    %testing plots 
    figure;
    scatter(wfPC1,wfPC2)
    xlabel("wfPC1")
    ylabel("wfPC2")

    figure;
    scatter (wfPC1,DS2_amplitudes)
    xlabel("wfPC1")
    ylabel("DS2 Amplitudes")

    figure;
    scatter (wfPC2,DS2_amplitudes)
    xlabel("wfPC2")
    ylabel("DS2 Amplitudes")

    figure;
    scatter (rateChange,DS2_amplitudes)
    xlabel("firing rate change")
    ylabel("DS2 Amplitudes")

    figure;
    scatter (burstIndex,DS2_amplitudes)
    xlabel("burstIndex")
    ylabel("DS2 Amplitudes")
    
    figure;
    scatter3(rateChange,DS2_amplitudes, wfPC1)
    xlabel("firing rate change")
    ylabel("DS2 Amplitudes")
    zlabel("wfPC1")

    figure()
    histogram(wfPC1)
    xlabel('w-PC1')
    ylabel('Cell Count')
    title('Histogram of PC1')
    
    figure()
    histogram(wfPC2)
    xlabel('w-PC2')
    ylabel('Cell Count')
    title('Histogram of PC2')


end 

function [tetShankChan] = makeTetShankChan(spatData,max_wf_chan)
    %make tetrode labels from CellIDs in spatData - can probably write
        %something else for this 
        cellInfo = getCellInfo(spatData);
        tet = cellInfo(:,2);
        % create channel inex per tetrode 'tet_index' which is a tetrode by
        % channel identity array
        tet_index = zeros(8, 4); % create an 8x4 matrix of zeros
        tet_index(1:2:end, :) = repmat(1:8:25, 4, 1).' + repmat(0:2:6, 4, 1); % fill odd-numbered rows with odd numbers
        tet_index(2:2:end, :) = repmat(2:8:26, 4, 1).' + repmat(0:2:6, 4, 1); % fill even-numbered rows with even numbers
        overlap = repmat (0:8:24,4,1).' + repmat (5:8,4,1);%overlap tetrodes are the bottom 4 contacts which are closer to each other 
        tet_index = [tet_index; overlap]; % + 1 overlap per octrode
        
        tetShankChan = zeros(length(tet),3);
        % this needs to be per probe type
        for it_tet = 1: length(tet)
            tetShankChan(it_tet,1) = tet(it_tet);
            tetShankChan(it_tet,3) = tet_index(tet(it_tet),max_wf_chan(it_tet));
            if tet(it_tet) == 1 || tet(it_tet) == 2 || tet(it_tet) == 9
                tetShankChan(it_tet,2) = 1;
            elseif tet(it_tet) == 3 || tet(it_tet) == 4 || tet(it_tet) == 10
                tetShankChan(it_tet,2) = 2;
            elseif tet(it_tet) == 5 || tet(it_tet) == 6 || tet(it_tet) == 11
                tetShankChan(it_tet,2) = 3;
            elseif tet(it_tet) == 7 || tet(it_tet) == 8 || tet(it_tet) == 12
                tetShankChan(it_tet,2) = 4;
            end
        end
end 

function [wfPC1,wfPC2] = waveformPCA(DG_ExCluster,waveforms)

    ex_waveforms = waveforms(DG_ExCluster);
    pca_data = [];
    for itEx = 1: length(ex_waveforms)
         if ~isnan(cell2mat(ex_waveforms(itEx))) % && length(cell2mat(ex_waveforms(itEx))) == 97 %changed from 50 don't know why this is here 
             pca_data = [pca_data; interp1(1:97, cell2mat(ex_waveforms(itEx)),1:0.48:97,'spline')];%[pca_data; interp1(1:50, cell2mat(ex_waveforms(itEx)),1:0.48:50,'spline')]
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

    % second derrivative 
    diff_pca_data = [];
    for itPCA = 1: size(pca_data,1)
        diff_pca_data = [diff_pca_data; diff(pca_data(itPCA,20:100),2)]; %inserted 1:80 for the .8 ms time window
    end
    [coeff,score] =  pca(diff_pca_data);
    w_PC1 = coeff(:,1);
    w_PC2 = coeff(:,2);
    
    %try just plotting sections of the second derrivative to see if you get
    %a bimodality 

    % weight of the PCA
    wfPC1 = score(:,1);
    wfPC2 = score(:,2);

    %dot product of 2nd derrivative of upsampled waveforms with w-PC1 and
    %w-PC2.
%     wfPC1 = [];
%     wfPC2   = [];
%  
%     for itPD = 1: size(diff_pca_data,1)
%          wfPC1 = [wfPC1; dot(diff_pca_data(itPD,:).',w_PC1)];
%          wfPC2 = [wfPC2; dot(diff_pca_data(itPD,:).',w_PC2)];
%     end
end