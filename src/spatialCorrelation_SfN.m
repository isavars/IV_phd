function spatialCorrelation_SfN(spatData)
%SPATIAL CORRELATION the goal is to have a remapping outcome measure like
%rate overlap but using spatial overlap. 
%   1.uses map_spatialcorr and rates_transform to make this score 
%   2.organized using age binning and indexes for environment.


    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env; 
    rMap = spatData.rMap;
    rMap1stHalf = spatData.rMap1stHalf;
    rMap2ndHalf = spatData.rMap2ndHalf;
    SI_spat = spatData.SI_spat;
    SpkTs = spatData.SpkTs;
    
    meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
%     SpkTs(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    
%gather age data from cellInfo 

    cellInfo = getCellInfo(spatData);
    age = cellInfo(:,3); 
    rat = cellInfo(:,1); 
    corecorded = cellInfo(:,4);
    
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

%make clusters
[cluster1,cluster2,cluster3,cluster4,cluster5] = makeClusters(spatData);

load('PCA2_clusters_from_WF1_AMR_BI_CCC.mat', 'PCA2_clusters', 'DG_ExCluster') %,'CA3_ExCluster') 

    %make clusters from PCA2_clusters keeping old naming convention for
    %convenience in running the old code. 
    cluster2 =[];
    cluster3 =[]; %keeping og naming convention for a second to see if this code runs 
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 1
            cluster3 = [cluster3;DG_ExCluster(ii)]; %for WF1_AMR_BI_CCC - 1 is mc and 2 is gc 
        elseif PCA2_clusters(ii) == 2
            cluster2 = [cluster2;DG_ExCluster(ii)];% for WF1_AMR_BI_DS2 - 1 is gc and 2 is mc 

        end
    end 
      
% % animal binning 
%    ratBins = [804 889];
%    for itRat = 1:size(ratBins,1)
%        indRat = rat>=ratBins(itRat,1) & rat<=ratBins(itRat,2);
%        meanRate = meanRate(indRat,:); 
%        burstIndex = burstIndex(indRat,:);

Age =[]; 
Environment = [];
spatCorrs = [];
FAMfieldNumC2 = [];
FAMfieldNumC3 = [];

NOVfieldNumC2 = [];
NOVfieldNumC3 = [];

%   age binning 

    ageBins   =  [17 20]; %list of age bins each spanning from col1:col2
    
    for itAge=1:size(ageBins,1)
        
        indAge = age>=ageBins(itAge,1) & age<=ageBins(itAge,2); % index for current age bin
        %meanRate = meanRate(indAge,:); 
        %burstIndex = burstIndex(indAge,:);

        %get indices for new clusters
        ages_indexes_in_spatData = [];
        for it_indAge =1: length(indAge) 
            if true(indAge(it_indAge))
                ages_indexes_in_spatData = [ages_indexes_in_spatData;it_indAge];
            else
            end
        end

       Newcluster2 = ismember(cluster2,ages_indexes_in_spatData);
       cluster2 = cluster2(Newcluster2);
       Newcluster3 = ismember(cluster3,ages_indexes_in_spatData);
       cluster3 = cluster3(Newcluster3);
    
        % make clusters (i want this to be its own function)
        meanMeanRate = [];
        meanBurstIndex = [];
        for it_gm = 1: length (meanRate)
             if SI_spat(it_gm,FamInd(it_gm,1))  >= 0.3 || SI_spat(it_gm,FamInd(it_gm,2))  >= 0.3 || SI_spat(it_gm,NovInd(it_gm))  >= 0.3 %
                meanMeanRate = [meanMeanRate; meanRate(it_gm,5)];
                meanBurstIndex = [meanBurstIndex; burstIndex(it_gm,5)];
             end
        end
        

        CellCount1 = length(cluster1);
        CellCount2 = length(cluster2);
        CellCount3 = length(cluster3);
        CellCount4 = length(cluster4);
        
        % make squares from circles 

        for itCr = 1: length (rMap)
            rMap(itCr,NovInd(itCr));
            Circle = rMap(itCr,NovInd(itCr));
            Square = cell2mat(rMap(itCr,FamInd(itCr,1)));
            NewSquare(itCr,:) = rates_transform(Circle,1,Square,8);
        end
        
        %get r values for all comparisons 

            FAMSpatCorrC1 = [];        NOVSpatCorrC1 = [];
            FAMSpatCorrC2 = [];        NOVSpatCorrC2 = [];
            FAMSpatCorrC3 = [];        NOVSpatCorrC3 = [];
            FAMSpatCorrC4 = [];        NOVSpatCorrC4 = [];
            SELFSpatCorrC1 = [];        
            SELFSpatCorrC2 = [];       
            SELFSpatCorrC3 = [];       
            SELFSpatCorrC4 = [];        
        
            for itRM = 1: length (cluster1)
                FAMSpatCorrC1(itRM, :) = map_spatialcorr(cell2mat(rMap(cluster1(itRM),FamInd(cluster1(itRM),1))), cell2mat(rMap(cluster1(itRM),FamInd(cluster1(itRM),2))));
                NOVSpatCorrC1(itRM, :) = map_spatialcorr(cell2mat(rMap(cluster1(itRM),FamInd(cluster1(itRM),2))), cell2mat(NewSquare(cluster1(itRM),1)));
                SELFSpatCorrC1(itRM, :) = map_spatialcorr(cell2mat(rMap1stHalf(cluster1(itRM),FamInd(cluster1(itRM),1))), cell2mat(rMap2ndHalf(cluster1(itRM),FamInd(cluster1(itRM),1))));
            end 
            for itRM = 1: length (cluster2)  
                FAMSpatCorrC2(itRM, :) = map_spatialcorr(cell2mat(rMap(cluster2(itRM),FamInd(cluster2(itRM),1))), cell2mat(rMap(cluster2(itRM),FamInd(cluster2(itRM),2))));
                NOVSpatCorrC2(itRM, :) = map_spatialcorr(cell2mat(rMap(cluster2(itRM),FamInd(cluster2(itRM),2))), cell2mat(NewSquare(cluster2(itRM),1)));
                SELFSpatCorrC2(itRM, :) = map_spatialcorr(cell2mat(rMap1stHalf(cluster2(itRM),FamInd(cluster2(itRM),1))), cell2mat(rMap2ndHalf(cluster2(itRM),FamInd(cluster2(itRM),1))));
            end
            for itRM = 1: length (cluster3) 
                FAMSpatCorrC3(itRM, :) = map_spatialcorr(cell2mat(rMap(cluster3(itRM),FamInd(cluster3(itRM),1))), cell2mat(rMap(cluster3(itRM),FamInd(cluster3(itRM),2))));
                NOVSpatCorrC3(itRM, :) = map_spatialcorr(cell2mat(rMap(cluster3(itRM),FamInd(cluster3(itRM),2))), cell2mat(NewSquare(cluster3(itRM),1)));
                SELFSpatCorrC3(itRM, :) = map_spatialcorr(cell2mat(rMap1stHalf(cluster3(itRM),FamInd(cluster3(itRM),1))), cell2mat(rMap2ndHalf(cluster3(itRM),FamInd(cluster3(itRM),1))));
            end
            for itRM = 1: length (cluster4)  
                FAMSpatCorrC4(itRM, :) = map_spatialcorr(cell2mat(rMap(cluster4(itRM),FamInd(cluster4(itRM),1))), cell2mat(rMap(cluster4(itRM),FamInd(cluster4(itRM),2))));
                NOVSpatCorrC4(itRM, :) = map_spatialcorr(cell2mat(rMap(cluster4(itRM),FamInd(cluster4(itRM),2))), cell2mat(NewSquare(cluster4(itRM),1)));
                SELFSpatCorrC4(itRM, :) = map_spatialcorr(cell2mat(rMap1stHalf(cluster4(itRM),FamInd(cluster4(itRM),1))), cell2mat(rMap2ndHalf(cluster4(itRM),FamInd(cluster4(itRM),1))));
            end
            

    %make means 
    
    FAMOverlapC1 = mean(rmmissing(FAMSpatCorrC1));
    NOVOverlapC1 = mean(rmmissing(NOVSpatCorrC1));
    FAMOverlapC2 = mean(rmmissing(FAMSpatCorrC2));
    NOVOverlapC2 = mean(rmmissing(NOVSpatCorrC2));
    FAMOverlapC3 = mean(rmmissing(FAMSpatCorrC3));
    NOVOverlapC3 = mean(rmmissing(NOVSpatCorrC3));
    FAMOverlapC4 = mean(rmmissing(FAMSpatCorrC4));
    NOVOverlapC4 = mean(rmmissing(NOVSpatCorrC4));
    
    SELFOverlapC1 = mean(rmmissing(SELFSpatCorrC1));
    SELFOverlapC2 = mean(rmmissing(SELFSpatCorrC2));
    SELFOverlapC3 = mean(rmmissing(SELFSpatCorrC3));
    SELFOverlapC4 = mean(rmmissing(SELFSpatCorrC4));

    
    % make errors
    FAMErrC1 = std(rmmissing(FAMSpatCorrC1))/sqrt(length(cluster1)); NOVErrC1 = std(rmmissing(FAMSpatCorrC1))/sqrt(length(cluster1));    
    FAMErrC2 = std(rmmissing(FAMSpatCorrC2))/sqrt(length(cluster2)); NOVErrC2 = std(rmmissing(FAMSpatCorrC2))/sqrt(length(cluster2));
    FAMErrC3 = std(rmmissing(FAMSpatCorrC3))/sqrt(length(cluster3)); NOVErrC3 = std(rmmissing(FAMSpatCorrC3))/sqrt(length(cluster3));    
    FAMErrC4 = std(rmmissing(FAMSpatCorrC4))/sqrt(length(cluster4)); NOVErrC4 = std(rmmissing(FAMSpatCorrC4))/sqrt(length(cluster4)); 
    
    %make figures             
        figure ()
        y = [ FAMOverlapC2 NOVOverlapC2; FAMOverlapC3 NOVOverlapC3];%;FAMOverlapC1 NOVOverlapC1; FAMOverlapC4 NOVOverlapC4];
        xticks = ({strcat('GCs : ',num2str(CellCount2)),strcat('MCs : ',num2str(CellCount3))});%,strcat('INs : ',num2str(CellCount1)),strcat('C4: ',num2str(CellCount4))});
        errors = [FAMErrC2 NOVErrC2; FAMErrC3 NOVErrC3];%; FAMErrC1  NOVErrC1;FAMErrC4 NOVErrC4];
        errorBars = gra_groupedbars(y, errors);
        errorBars = set(gca,'xticklabels', xticks);
        title(strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2))));
                ylabel('Spatial Correlation');
                xlabel('Cell Type: N');
        
%                 ylim([-0.1 0.35]);
                
%         figure ()
%         x = categorical({'C1','C2','C3','C4'});
%         y = [SELFOverlapC1 ; SELFOverlapC2; SELFOverlapC3; SELFOverlapC4];
%         SpatialCorr= bar(x,y);
%             SpatialCorr(1).FaceColor = 'flat';
%             SpatialCorr(1).CData(1,:) = [1 0 0];
%             SpatialCorr(1).CData(2,:) = [0 0 1];
%             SpatialCorr(1).CData(3,:) = [0 1 0];
%             SpatialCorr(1).CData(4,:) = [0 0 0];    
%                 title(strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2))));
%                 ylabel('Intratrial Spatial Correlation');
%                 ylim([-0.1 0.35]);
    
%     %make parent figure with axes to be filled in
%     
%         hFig = gra_multiplot(1, 4, 'figborder', [2 1 1 1]);
%         axArr = getappdata(hFig, 'axesHandles' ); % makes the axes 
% 
%            g1 = categorical({'C1F','C1N'});
%            g2 = categorical({'C2F','C2N'});
%            g3 = categorical({'C3F','C3N'});
%            g4 = categorical({'C4F','C4N'});           
%            boxplot(axArr(1,1),[FAMSpatCorrC1, NOVSpatCorrC1], g1)
%            boxplot(axArr(1,2),[FAMSpatCorrC2, NOVSpatCorrC2], g2) 
%            boxplot(axArr(1,3),[FAMSpatCorrC3, NOVSpatCorrC3], g3)
%            boxplot(axArr(1,4),[FAMSpatCorrC4, NOVSpatCorrC4], g4)
%            ylabel(axArr(1,1), 'Spatial Correlation')
%             title(axArr(1,1), strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2))));
         meanRate = spatData.meanRate;
         burstIndex = spatData.burstIndex;
         
         FAMSpatCorrC2 = atanh(FAMSpatCorrC2);% to make it parametric for mixed anova
         NOVSpatCorrC2 = atanh(NOVSpatCorrC2);
         FAMSpatCorrC3 = atanh(FAMSpatCorrC3);% to make it parametric for mixed anova
         NOVSpatCorrC3 = atanh(NOVSpatCorrC3);
         % stats
        spatCorrs = [spatCorrs; FAMSpatCorrC2 ;NOVSpatCorrC2];%[spatCorrs; NOVSpatCorrC2 -  FAMSpatCorrC2];%[rmmissing(FAMSpatCorrC3) ;rmmissing(NOVSpatCorrC3)];
        Age = [Age;repmat(ageBins(itAge,1),length(FAMSpatCorrC2),1);repmat(ageBins(itAge,1),length(NOVSpatCorrC2),1)];
        Environment = [Environment; repmat('Fam',length(FAMSpatCorrC2),1);repmat('Nov', length(NOVSpatCorrC2),1)];
        
        % make number of fields (spatial only) %also make for fam2 and nov can get visual
        % of remapping 

%             fieldNumC2 = [];
            fieldCount = 0;
            for itFC = 1: length (cluster2)
                map_rate_thr = cell2mat(rMap(cluster2(itFC),FamInd(itFC,1))) >= meanRate(cluster2(itFC),FamInd(itFC,1));
                stats = regionprops(map_rate_thr, 'Area');
                Areas = cell2mat(struct2cell(stats));
                for itA = 1:length(Areas)
                    if Areas(itA) >= 50
                        fieldCount = fieldCount + 1;
                    end
                end
                FAMfieldNumC2 =[FAMfieldNumC2; fieldCount];
                fieldCount = 0;
            end
            fieldCount = 0;
            for itFC = 1: length (cluster2)
                map_rate_thr = cell2mat(rMap(cluster2(itFC),NovInd(itFC))) >= meanRate(cluster2(itFC),NovInd(itFC));
                stats = regionprops(map_rate_thr, 'Area');
                Areas = cell2mat(struct2cell(stats));
                for itA = 1:length(Areas)
                    if Areas(itA) >= 20
                        fieldCount = fieldCount + 1;
                    end
                end
                NOVfieldNumC2 =[NOVfieldNumC2; fieldCount];
                fieldCount = 0;
            end
            
%             fieldNumC3 = [];
            fieldCount = 0;
            for itFC = 1: length (cluster3)
                map_rate_thr = cell2mat(rMap(cluster3(itFC),FamInd(itFC,1))) >= meanRate(cluster3(itFC),FamInd(itFC,1));
                stats = regionprops(map_rate_thr, 'Area');
                Areas = cell2mat(struct2cell(stats));
                for itA = 1:length(Areas)
                    if Areas(itA) >= 16
                        fieldCount = fieldCount + 1;
                    end
                end
                FAMfieldNumC3 =[FAMfieldNumC3; fieldCount];
                fieldCount = 0;
            end
            fieldCount = 0;
            for itFC = 1: length (cluster3)
                map_rate_thr = cell2mat(rMap(cluster3(itFC),NovInd(itFC))) >= meanRate(cluster3(itFC),NovInd(itFC));
                stats = regionprops(map_rate_thr, 'Area');
                Areas = cell2mat(struct2cell(stats));
                for itA = 1:length(Areas)
                    if Areas(itA) >= 20
                        fieldCount = fieldCount + 1;
                    end
                end
                NOVfieldNumC3 =[NOVfieldNumC3; fieldCount];
                fieldCount = 0;
            end
    end
   %end             
%     [p,tbl,stats]= kruskalwallis(spatCorrs, Age);
%     [p,tbl,stats]=anovan(spatCorrs,{Age,Environment},'model','interaction', 'varnames',{'Age','Environment'})
    figure()
    histogram(FAMfieldNumC2)
    figure()
    histogram(FAMfieldNumC3)
    
    h= kstest2(FAMfieldNumC2,FAMfieldNumC3);
    h= kstest2(NOVfieldNumC2,NOVfieldNumC3);
%     gc_spatRem = [FAMfieldNumC2,NOVfieldNumC2];
%     mc_spatRem = [FAMfieldNumC3,NOVfieldNumC3];
%     save('spatRemData.mat','gc_spatRem','mc_spatRem')
end


