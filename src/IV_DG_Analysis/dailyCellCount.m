function dailyCellCount()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load ('allRats_DGCA3_spatData.mat', 'spatData')

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env; 
    SI_spat = spatData.SI_spat;
    wf_means = spatData.wf_means;
    nSpks = spatData.nSpks;
    rMap = spatData.rMap;
    
    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    
%gather age data from cellInfo 

    cellInfo = getCellInfo();
    corecorded = cellInfo(:,4);
    age = cellInfo(:,3); 
    rat = cellInfo(:,1); 
    
%make indexes for environment to use in comparisons with any file 
        maxFam = 2;
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
        awakeMeanRate(it_gma,:) = mean(meanRate(it_gma,(FamInd(1):FamInd(2))));
    end

    sleepMeanRate = meanRate (:,5);
    
    rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));
    
% values for table 

    Age = [];
    Rats = [];
    cellNum = [];
    fieldNum = [];
    
%   age binning 

    ageBins   = [18 19;20 21;22 23];  % list of age bins each spanning from col1:col2
    
    for itAge=1:size(ageBins,1)
        
        indAge = age>=ageBins(itAge,1) & age<=ageBins(itAge,2); % index for current age bin
        meanRate = meanRate(indAge,:); 
        burstIndex = burstIndex(indAge,:);
        
        % make clusters (i want this to be its own function)
        meanMeanRate = [];
        meanBurstIndex = [];
        for it_gm = 1: length (meanRate)
                meanMeanRate = [meanMeanRate; meanRate(it_gm,5)];% will change when I sort out final clustering method 
                meanBurstIndex = [meanBurstIndex; burstIndex(it_gm,5)];
        end
        
        cluster1 = [];
        cluster2 = [];
        cluster3 = [];
        cluster4 = [];

        for jj = 1: length(meanMeanRate)
            if  (meanMeanRate(jj) >= 1) && (meanMeanRate(jj) <= 10) && (meanBurstIndex(jj) <= 0.05) %interneurons
                cluster1 = [cluster1; jj];
            elseif  (rateChange(jj) <= 0.3) && (meanMeanRate(jj) >= 0.01) && (meanMeanRate(jj) < 1.5)  && (meanBurstIndex(jj) >= 0.01) && (meanBurstIndex(jj) <= 0.15) %potential granule cells
                cluster2 = [cluster2; jj];
            elseif (meanMeanRate(jj) >= 0.02) && (meanBurstIndex(jj) >= 0.15) %potential mossy cells
                cluster3 = [cluster3; jj];
            else
                cluster4 = [cluster4; jj];
            end    
        end 
 
         % make # of cells in a cluster per rat per day 
        
        CellCount1 = length(cluster1);
        CellCount2 = length(cluster2);
        CellCount3 = length(cluster3);
        CellCount4 = length(cluster4);
        
        rats = rat(indAge,:);
        rats = rats(cluster2);
        
        P = repmat(ageBins(itAge,1), CellCount2,1);   
        
        indRat = diff(rats) ~= 0;  indRat = find(indRat == 1);
        if ~isempty(indRat)
            if length(indRat) >= 2 
                Num = repmat(indRat(1), indRat(1),1);
                Num = [Num;repmat(indRat(2) - indRat(1), indRat(2) - indRat(1),1); repmat(length(rats)- indRat(2),length(rats)- indRat(2),1)]; 
            else
                Num = [repmat(indRat, indRat,1); repmat(length(rats)- indRat,length(rats)- indRat,1)];
            end
        elseif isempty(indRat) 
            Num = repmat(CellCount2, CellCount2,1);
        end
        
% make number of fields 

    fieldNumC2 = [];
    fieldCount = 0;
    for itFC = 1: length (cluster2)
        map_rate_thr = cell2mat(rMap(cluster2(itFC),FamInd(1))) >= meanRate(cluster2(itFC),FamInd(1));
        stats = regionprops(map_rate_thr, 'Area');
        Areas = cell2mat(struct2cell(stats));
        for itA = 1:length(Areas)
            if Areas(itA) >= 20
                fieldCount = fieldCount + 1;
            end
        end
        fieldNumC2 =[fieldNumC2; fieldCount];
        fieldCount = 0;
    end
    
    fieldNumC3 = [];
    fieldCount = 0;
    for itFC = 1: length (cluster3)
        map_rate_thr = cell2mat(rMap(cluster3(itFC),FamInd(1))) >= meanRate(cluster3(itFC),FamInd(1));
        stats = regionprops(map_rate_thr, 'Area');
        Areas = cell2mat(struct2cell(stats));
        for itA = 1:length(Areas)
            if Areas(itA) >= 20
                fieldCount = fieldCount + 1;
            end
        end
        fieldNumC3 =[fieldNumC3; fieldCount];
        fieldCount = 0;
    end

    %interested in using pie charts to display info function is called pie 

        %cluster = [cluster1, cluster2, cluster3, cluster4];
%         figure()
%         boxplot(corecorded(cluster2));
%         title({'Cluster 2: ', num2str(ageBins(itAge,1)), '(', num2str(CellCount2), ')'})        
%         figure()
%         boxplot(corecorded(cluster3));
%         title({'Cluster 3: ', num2str(ageBins(itAge,1)), '(', num2str(CellCount3), ')'})
%         figure()
%         boxplot(corecorded(cluster4));
%         title({'Cluster 4: ', num2str(ageBins(itAge,1)), '(', num2str(CellCount4), ')'})
%         
%         corecordedMeans = [mean(corecorded(cluster2)), mean(corecorded(cluster3)),mean(corecorded(cluster4))];
%         bar(corecordedMeans);
        
        Rats = [Rats;rats];
        Age = [Age;P];
        cellNum = [cellNum;Num]; 
        rat = cellInfo(:,1); 
        meanRate = spatData.meanRate;
        burstIndex = spatData.burstIndex;
        
    end   

    %make table     
        varNames = {'Age', 'Rat', 'Granule Cell Number'}; 
        cellType = {'Interneuron'};
        dailyCellCount = table(Age,Rats,cellNum, 'VariableNames', varNames );
        unique(dailyCellCount)
        

        
        figure()
        c2 =boxplot(burstIndex(cluster2));
            ylabel('Burst Index')
            xlabel('Putative Granule Cells')
%         figure()
%         c3 =boxplot(burstIndex(cluster3));
%             ylabel('Burst Index')
%             xlabel('Putative Mossy Cells')
%         figure()
%         c2 =boxplot(meanRate(cluster2));
%             ylabel('Mean Rate(Hz)')
%             xlabel('Putative Granule Cells')
%             set(gca, 'yscale', 'log')
%         figure()
%         c3 =boxplot(meanRate(cluster3));
%             ylabel('Mean Rate(Hz)')
%             xlabel('Putative Mossy Cells')
%             set(gca, 'yscale', 'log')
end        
        


