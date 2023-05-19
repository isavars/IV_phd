function spatialCorrelation_mossy_SfN(spatData)
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
    
%     meanRate(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
%     burstIndex(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
%     SI_spat(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
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
        DiffInd = [];
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
                elseif strcmp(cast(spatData.env(itCell,itTrial),'char'),'diff')
                    DiffInd(itCell,itTrial) = itTrial;
                    DiffInd = nonzeros(DiffInd);
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

      
% % animal binning 
%    ratBins = [804 889];
%    for itRat = 1:size(ratBins,1)
%        indRat = rat>=ratBins(itRat,1) & rat<=ratBins(itRat,2);
%        meanRate = meanRate(indRat,:); 
%        burstIndex = burstIndex(indRat,:);

Age =[]; 
Environment = [];
spatCorrs = [];
SI_score =[];
Ages =[];
FAMfieldNumC2 = [];
FAMfieldNumC3 = [];

NOVfieldNumC2 = [];
NOVfieldNumC3 = [];

DIFFfieldNumC3 =[];

%   age binning 

    ageBins   =  [17 20; 22 31];  %[18 23]; list of age bins each spanning from col1:col2
    original_cluster = cluster3; 
    
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
    
%         % make clusters (i want this to be its own function)
%         meanMeanRate = [];
%         meanBurstIndex = [];
%         for it_gm = 1: length (meanRate)
%              if SI_spat(it_gm,FamInd(it_gm,1))  >= 0.3 || SI_spat(it_gm,FamInd(it_gm,2))  >= 0.3 || SI_spat(it_gm,NovInd(it_gm))  >= 0.3 %
%                 meanMeanRate = [meanMeanRate; meanRate(it_gm,5)];
%                 meanBurstIndex = [meanBurstIndex; burstIndex(it_gm,5)];
%              end
%         end


        CellCount2 = length(cluster2);
        CellCount3 = length(cluster3);

        
        % make squares from circles 

        for itCr = 1: length (rMap)
            rMap(itCr,NovInd(itCr));
            Circle = rMap(itCr,NovInd(itCr));
            Square = cell2mat(rMap(itCr,FamInd(itCr,1)));
            NewSquare(itCr,:) = rates_transform(Circle,1,Square,8);
        end
        
        %get r values for all comparisons 


        FAMSpatCorrC3 = [];        
        NOVSpatCorrC3 = [];
        DIFFSpatCorrC3 =[];       
        SELFSpatCorrC3 = [];       
  
    
        for itRM = 1: length (cluster3) 
            FAMSpatCorrC3(itRM, :) = map_spatialcorr(cell2mat(rMap(cluster3(itRM),FamInd(cluster3(itRM),1))), cell2mat(rMap(cluster3(itRM),FamInd(cluster3(itRM),2))));
            NOVSpatCorrC3(itRM, :) = map_spatialcorr(cell2mat(rMap(cluster3(itRM),FamInd(cluster3(itRM),2))), cell2mat(NewSquare(cluster3(itRM),1)));
            SELFSpatCorrC3(itRM, :) = map_spatialcorr(cell2mat(rMap1stHalf(cluster3(itRM),FamInd(cluster3(itRM),1))), cell2mat(rMap2ndHalf(cluster3(itRM),FamInd(cluster3(itRM),1))));
            DIFFSpatCorrC3(itRM,:) = map_spatialcorr(cell2mat(rMap(cluster3(itRM),FamInd(cluster3(itRM),2))),cell2mat(rMap(cluster3(itRM),DiffInd(cluster3(itRM),1))));
        end

            

        %make means 
            
    
        FAMOverlapC3 = mean(rmmissing(FAMSpatCorrC3));
        NOVOverlapC3 = mean(rmmissing(NOVSpatCorrC3));    
        DIFFOverlapC3 = mean(rmmissing(DIFFSpatCorrC3)); 
        SELFOverlapC3 = mean(rmmissing(SELFSpatCorrC3));
      
        
        % make errors
         FAMErrC3 = std(rmmissing(FAMSpatCorrC3))/sqrt(length(cluster3)); NOVErrC3 = std(rmmissing(FAMSpatCorrC3))/sqrt(length(cluster3));    
         DIFFErrC3 = std(rmmissing(DIFFSpatCorrC3))/sqrt(length(cluster3));
    
        %make figures             
        figure ()
        x = categorical({'FamXFam','FamXNov1','FamXNov2'});
        y = [FAMOverlapC3 DIFFOverlapC3 NOVOverlapC3];%;FAMOverlapC1 NOVOverlapC1; FAMOverlapC4 NOVOverlapC4];
        xticks = ({strcat('MCs : ',num2str(CellCount3))});%,strcat('INs : ',num2str(CellCount1)),strcat('C4: ',num2str(CellCount4))});
        errors = [FAMErrC3 DIFFErrC3 NOVErrC3];%; FAMErrC1  NOVErrC1;FAMErrC4 NOVErrC4];
%         errorBars = gra_groupedbars(y, errors);
        SpatialCorr= bar(x,y);
        hold on
            er = errorbar(x,y,errors,errors);    
            er.Color = [0 0 0];                            
            er.LineStyle = 'none';  
        hold off
%         errorBars = set(gca,'xticklabels', xticks);
        title(strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2))));
                ylabel('Spatial Correlation');
                xlabel('Cell Type: N');
                ylim([0 0.6]);
                
         meanRate = spatData.meanRate;
         burstIndex = spatData.burstIndex;
         cluster3= original_cluster;

%          FAMSpatCorrC3 = atanh(FAMSpatCorrC3);% to make it parametric for mixed anova
%          NOVSpatCorrC3 = atanh(NOVSpatCorrC3);
%          DIFFSpatCorrC3 = atanh(DIFFSpatCorrC3);
         
         meanSIs = SI_spat(cluster3);
         % stats
         SI_score = [SI_score;meanSIs];
         Ages = [Ages;repmat(ageBins(itAge,1),length(cluster3),1)];
         
        Age = [Age;itAge*ones(length(FAMSpatCorrC3)*3,1)];%;repmat(ageBins(itAge,1),length(cluster3),1);repmat(ageBins(itAge,1),length(cluster3),1)];% *2 removed for KW atempt 
        spatCorrs = [spatCorrs;FAMSpatCorrC3 ;DIFFSpatCorrC3; NOVSpatCorrC3];%[rateOverlaps; novOverlap2 - famOverlap2];%[famOverlap2 ;novOverlap2];%[rateOverlaps;famOverlap3 ,novOverlap3]%
        Environment = [Environment;ones(length( NOVSpatCorrC3),1); 2*ones(length( NOVSpatCorrC3),1);3*ones(length( NOVSpatCorrC3),1)];
         %spatCorrs = [spatCorrs;FAMSpatCorrC3; DIFFSpatCorrC3; NOVSpatCorrC3];%[spatCorrs; NOVSpatCorrC2 -  FAMSpatCorrC2];%[rmmissing(FAMSpatCorrC3) ;rmmissing(NOVSpatCorrC3)];
         %Age = [Age;repmat(ageBins(itAge,1),length(spatCorrs),1)];%repmat(ageBins(itAge,1),length(cluster3),1);repmat(ageBins(itAge,1),length(cluster3),1)];
         %Environment = [Environment; repmat('Fam',length(FAMSpatCorrC3),1);repmat('Diff', length(DIFFSpatCorrC3),1);repmat('Nov', length(NOVSpatCorrC3),1)];
         CellCount3 = 0;
        
 end       
       
    %t-test for spatial score and age
       %[p,tbl] = anova1(SI_score,Ages);
       [p,tbl,stats]=anovan(spatCorrs,{Age,Environment},'model','interaction', 'varnames',{'Age','Environment'});
        c = multcompare(stats)
%     save('spatRemData.mat','gc_spatRem','mc_spatRem')

end 


