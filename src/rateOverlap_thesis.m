function rateOverlap_thesis(data, electrode_positions)
%RATEOVERLAP is based on the Leutgeb et al., 2007 analysis: 
%   1. Divide (for each cell) the mean firing rate in the less active
%   enclosure by the mean firing rate in the more active enclosure.
%   2. Average the rates of the active cells for the cell populaiton data.
%   The resulting overlap scores are 1 if the cell populaiton is compared
%   with itself and 0 if a completely difrerent cell population is active.
%   3. should probably do the same stats as in the Leutgeb paper or
%   laurenz's number nSpks in each environment to be considered.  

% load data as input select the .mat file with the data you want and run
% with spatData as the input. 

load(data, 'spatData')
load(electrode_positions, 'elePos') %necessary for now for running class_cells but might change in the future 

    meanRate = spatData.meanRate;
    burstIndex = spatData.burstIndex;
    env = spatData.env; 
    SI_spat = spatData.SI_spat;
    nSpks = spatData.nSpks;
    
    meanRate(~(env=='fam' | env=='nov' | env=='diff' | env=='sleep')) = NaN;
    burstIndex(~(env=='fam' | env=='nov'| env=='diff'  | env=='sleep')) = NaN;
    SI_spat(~(env=='fam' | env=='nov' | env=='diff' | env=='sleep')) = NaN;


%gather age data from cellInfo %figure out how to change this 

    cellInfo = getCellInfo(spatData);
    age = cellInfo(:,3);
    rat = cellInfo(:,1); 
    
%make indexes for environment to use in comparisons with any file 
    maxFam = 2;
    FamInd = nan(length(meanRate),maxFam); 
    FamIndT = [];     
    NovInd = [];
    DiffInd = [];
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
        SI_spat(it_gma,:) = max(SI_spat(it_gma,(1:6)));
    end

    sleepMeanRate = meanRate (:,6);
    
    rateChange = ((awakeMeanRate-sleepMeanRate)./(awakeMeanRate+sleepMeanRate));


%get clusters from class_cells (can also save final classification for this
%to not rely on clusters being made every time. 

    [DS2_amplitude,PCA2_clusters, DG_ExCluster, InCluster] = class_cells(data, electrode_positions); %check which outputs you need 
    
    
Age =[];
rateOverlaps =[];
Environment = [];
%   age binning 

    ageBins   =  [21 24; 25 31];  %list of age bins each spanning from col1:col2

%     clusters = {cluster2,cluster3};%array of cluster names 
%     original_cluster = cluster3;
    cluster2 =[];
    cluster3 =[]; %keeping og naming convention for a second to see if this code runs 
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 1 
            cluster3 = [cluster3;DG_ExCluster(ii)]; 
        elseif PCA2_clusters(ii) == 2
            cluster2 = [cluster2;DG_ExCluster(ii)];

        end
    end 

    original_cluster = cluster3;
    
    for itAge=1:size(ageBins,1)
        
        indAge = age>=ageBins(itAge,1) & age<=ageBins(itAge,2); % index for current age bin

        meanMeanRate = [];
        meanBurstIndex = [];
        for it_gm = 1: length (meanRate)
            if nSpks(it_gm,FamInd(it_gm,1))  >= 10 || nSpks(it_gm,FamInd(it_gm,2))  >= 10 || nSpks(it_gm,NovInd(it_gm))  >= 10
                meanMeanRate(it_gm) = mean(meanRate(it_gm,6));% indexing here needs to generalize to all datasets 
                meanBurstIndex(it_gm) = mean(burstIndex(it_gm,6));
            end
        end

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
      

% %          residual = novOverlap2 - NOVOverlapC2
%          figure;
%          histogram (residual)

%          famOverlap2 = atanh(famOverlap2);% to make it parametric for
%          mixed anova
%          novOverlap2 = atanh(novOverlap2);
         
%          
        famOverlap3 = [];
        novOverlap3 = [];
        diffOverlap3= [];
         for itC1 = 1: length(cluster3)
%              if (nSpks(cluster3(itC1),1) || nSpks(cluster3(itC1),2) || nSpks(cluster3(itC1),3)) <= 1
%              end
             if meanRate(cluster3(itC1),FamInd(cluster3(itC1),1)) > meanRate(cluster3(itC1),FamInd(cluster3(itC1),2))
                famOverlap3(itC1) = meanRate (cluster3(itC1),FamInd(cluster3(itC1),2)) / meanRate(cluster3(itC1),FamInd(cluster3(itC1),1));
             elseif meanRate(cluster3(itC1),FamInd(cluster3(itC1),1)) < meanRate(cluster3(itC1),FamInd(cluster3(itC1),2))
                famOverlap3(itC1) = meanRate (cluster3(itC1),FamInd(cluster3(itC1),1)) / meanRate(cluster3(itC1),FamInd(cluster3(itC1),2)); %this isn't going to work with all data some have nan in position 1 so i have to use 2vs 4 
             end
             if meanRate(cluster3(itC1),FamInd(cluster3(itC1),2)) > meanRate(cluster3(itC1),NovInd(cluster3(itC1)))
                novOverlap3(itC1) = meanRate (cluster3(itC1),NovInd(cluster3(itC1))) / meanRate(cluster3(itC1),FamInd(cluster3(itC1),2));
             elseif meanRate(cluster3(itC1),FamInd(cluster3(itC1),2)) < meanRate(cluster3(itC1),NovInd(cluster3(itC1)))
                novOverlap3(itC1) = meanRate (cluster3(itC1),FamInd(cluster3(itC1),2)) / meanRate(cluster3(itC1),NovInd(cluster3(itC1)));
             end
             if meanRate(cluster3(itC1),FamInd(cluster3(itC1),2)) > meanRate(cluster3(itC1),DiffInd(cluster3(itC1)))
                diffOverlap3(itC1) = meanRate (cluster3(itC1),DiffInd(cluster3(itC1))) / meanRate(cluster3(itC1),FamInd(cluster3(itC1),2));
             elseif meanRate(cluster3(itC1),FamInd(cluster3(itC1),2)) < meanRate(cluster3(itC1),DiffInd(cluster3(itC1)))
                diffOverlap3(itC1) = meanRate (cluster3(itC1),FamInd(cluster3(itC1),2)) / meanRate(cluster3(itC1),DiffInd(cluster3(itC1)));
             end
         end
         famOverlap3 = rmmissing(transpose(famOverlap3));
         novOverlap3 = rmmissing(transpose(novOverlap3));
         diffOverlap3 = rmmissing(transpose(diffOverlap3));         
         FAMOverlapC3 = mean(famOverlap3);
         NOVOverlapC3 = mean(novOverlap3);
         DIFFOverlapC3 = mean(diffOverlap3);
         FAMErrC3 = std(famOverlap3)/sqrt(length(cluster3));
         NOVErrC3 = std(novOverlap3)/sqrt(length(cluster3));
         DIFFErrC3 = std(diffOverlap3)/sqrt(length(cluster3));
         
%          famOverlap3 = atanh(famOverlap3);
%          novOverlap3 = atanh(novOverlap3);
% % 
%          

        CellCount3 = length(cluster3);


        figure()
        x = categorical({'FamXFam','FamXNov1','FamXNov2'});
        y = [FAMOverlapC3 DIFFOverlapC3 NOVOverlapC3];%; FAMOverlapC1 NOVOverlapC1; FAMOverlapC4 NOVOverlapC4];
        errors = [ FAMErrC3 DIFFErrC3 NOVErrC3]; %; FAMErrC1  NOVErrC1;FAMErrC4 NOVErrC4];
%         xticks = ({strcat('GCs (',num2str(CellCount2),')'),strcat('MCs (',num2str(CellCount3),')')}); %strcat('INs (',num2str(CellCount1),')'),,strcat('C4(',num2str(CellCount4),')')});
%         errorBars = gra_groupedbars(y, errors);
%         errorBars = set(gca,'xticklabels', xticks);
        RateOver= bar(x,y);
        hold on
            er = errorbar(x,y,errors,errors);    
            er.Color = [0 0 0];                            
            er.LineStyle = 'none';  
        hold off
        title(strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2)))); %r',num2str(ratBins(itRat,1)),': put this back in when you fix rat binning 
            ylabel('Rate Overlap Score')
            xlabel('Cell Type (N)')
        
         figure()
            y = mean(SI_spat(cluster3))
            errorC3 = transpose(std(SI_spat(cluster3))/sqrt(length(cluster3)));
            errors = errorC3
        %         xticks = ({strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2)))}); 
        %         errorBars = gra_groupedbars(y, errors);
        %         errorBars = set(gca,'xticklabels', xticks);
            bar(y)
            title(strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2)))); 
                ylabel('Spatial Information Score')
                xlabel('Age')

        Age = [Age;itAge*ones(length(famOverlap3)*3,1)];%;repmat(ageBins(itAge,1),length(cluster3),1);repmat(ageBins(itAge,1),length(cluster3),1)];% *2 removed for KW atempt 
        rateOverlaps = [rateOverlaps;famOverlap3 ;diffOverlap3;novOverlap3];%[rateOverlaps; novOverlap2 - famOverlap2];%[famOverlap2 ;novOverlap2];%[rateOverlaps;famOverlap3 ,novOverlap3]%
        Environment = [Environment;ones(length(famOverlap3),1); 2*ones(length(famOverlap3),1);3*ones(length(famOverlap3),1)];
         %Environment = [repmat('Fam',length(famOverlap3),1);repmat('diff',length(cluster3),1);repmat('Nov', length(novOverlap3),1)];
%         Environment = [Environment;repmat('Fam',length(cluster3),1);repmat('diff',length(cluster3),1);repmat('Nov', length(cluster3),1)];
        %Clusters = {'C2';'C2';'C3';'C3';'C4';'C4'};


        meanRate = spatData.meanRate;
        burstIndex = spatData.burstIndex;
        cluster3 = original_cluster;

    end  
   %end
%   [repmat('Fam',length(cluster3),1);repmat('diff',length(cluster3),1);repmat('Nov', length(cluster3),1)]

  
    [p,tbl,stats]=anovan(rateOverlaps,{Age,Environment},'model','interaction', 'varnames',{'Age','Environment'});
   
%     [p,tbl,stats]= kruskalwallis(rateOverlaps, Age);
    c = multcompare(stats)


    

end        

