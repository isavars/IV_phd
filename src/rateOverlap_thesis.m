function rateOverlap_thesis(data, cell_clusters)
%RATEOVERLAP is based on the Leutgeb et al., 2007 analysis: 
%   1. Divide (for each cell) the mean firing rate in the less active
%   enclosure by the mean firing rate in the more active enclosure.
%   2. Average the rates of the active cells for the cell population data.
%   The resulting overlap scores are 1 if the cell populaiton is compared
%   with itself and 0 if a completely difrerent cell population is active.
%   3. should probably do the same stats as in the Leutgeb paper or
%   laurenz's number nSpks in each environment to be considered.  

% INPUT - spatData from the dataset you want to test and the clusters from
% saved class_cells output. 

%load spatial data for all cells and gather necessary parameters 
load(data, 'spatData')

    meanRate = spatData.meanRate;
    env = spatData.env; 
    SI_spat = spatData.SI_spat;
    nSpks = spatData.nSpks;
    
%produce age data from cellInfo 
    cellInfo = getCellInfo(spatData);
    age = cellInfo(:,3);

    
%make indexes for environment to use in comparisons with any file - this
%should probably be a function 'environment_inds'- and I want it to find
%the best famBox trials out of the three instead of the first 2 
        maxFam = 2;
        FamInd = nan(height(spatData),maxFam); 
        FamIndT = [];     
        NovInd = [];
        DiffInd = [];
        famCount = 0;
        
        for itCell= 1: height(spatData)
            for itTrial = 1: width(env)
                if contains(cast(env(itCell,itTrial),'char'),'fam')
                    FamIndT(itCell,itTrial) = itTrial;
                    famCount = famCount + 1;
                    if famCount <= maxFam
                         FamInd(itCell,famCount)= transpose(nonzeros(FamIndT(itCell,itTrial)));
                    end
                elseif strcmp(cast(env(itCell,itTrial),'char'),'nov')
                    NovInd(itCell,itTrial) = itTrial;
                    NovInd = nonzeros(NovInd);
                elseif strcmp(cast(env(itCell,itTrial),'char'),'diff')
                    DiffInd(itCell,itTrial) = itTrial;
                    DiffInd = nonzeros(DiffInd);
                end 
            end
            famCount = 0;
        end

      
%make spatial information score to be used for     
    for it_sp = 1: height (spatData)
        SI_spat(it_sp,:) = max(SI_spat(it_sp,:),[],'omitnan');
    end

%get clusters from class_cells output
load(cell_clusters, 'PCA2_clusters', 'DG_ExCluster','CA3_ExCluster')
     
    %make clusters from PCA2_clusters keeping old naming convention for
    %convenience in running the old code. 
    mossy =[];
    granule =[]; %keeping og naming convention for a second to see if this code runs 
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 1
            granule = [granule;DG_ExCluster(ii)]; %for WF1_AMR_BI_CCC_latest - 2 is gc and 1 is mc for WF1_AMR_BI_CCC_new - 2 is gc and 1 is mc 
        elseif PCA2_clusters(ii) == 2
            mossy = [mossy;DG_ExCluster(ii)];% for WF1_AMR_BI_DS2_new - 2 is gc and 1 is mc 

        end
    end 

    cluster3 = granule ;%granule;%CA3_ExCluster;

    %make agebins and loop through to get cluster data for each age bin
    Age =[];
    rateOverlaps =[];
    Environment = [];
    ageBins   =  [17 20; 21 31];  %list of age bins each spanning from col1:col2
   
    original_cluster = cluster3; %no idea why this is here and down below
    length(cluster3)
    for itAge=1:size(ageBins,1)
        
        indAge = age>=ageBins(itAge,1) & age<=ageBins(itAge,2); % index for current age bin
        
        %make new clusters which are only for the current age bin
        ages_indexes_in_spatData = [];
        for it_indAge =1: length(indAge) 
            if true(indAge(it_indAge))
                ages_indexes_in_spatData = [ages_indexes_in_spatData;it_indAge];
            else
            end
        end
        Newcluster3 = ismember(cluster3,ages_indexes_in_spatData);
        cluster3 = cluster3(Newcluster3);

      
         %do rate overlap for comparisons between the three environments 
         famOverlap3 = nan(1,length(cluster3));
         novOverlap3 = nan(1,length(cluster3));
         diffOverlap3= nan(1,length(cluster3));
         for itC1 = 1: length(cluster3)
              if any(nSpks(cluster3(itC1),1:5) > 50) %exclude cells that didn't really fire in any of the boxes            
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
              else
              end
         end
        famOverlap3 = famOverlap3.';
        novOverlap3 = novOverlap3.';
        diffOverlap3 = diffOverlap3.';
        
        % Find indices where any of the arrays have a NaN
        nanIndices = isnan(famOverlap3) | isnan(novOverlap3) | isnan(diffOverlap3);        
        % Remove those indices from each array
        famOverlap3(nanIndices) = [];
        novOverlap3(nanIndices) = [];
        diffOverlap3(nanIndices) = [];  

        FAMOverlapC3 = mean(famOverlap3);
        NOVOverlapC3 = mean(novOverlap3);
        DIFFOverlapC3 = mean(diffOverlap3);
        FAMErrC3 = std(famOverlap3)/sqrt(length(cluster3));
        NOVErrC3 = std(novOverlap3)/sqrt(length(cluster3));
        DIFFErrC3 = std(diffOverlap3)/sqrt(length(cluster3));

        %make Agebin titles 

        if itAge == 1
            AgeBin = 'Pre-wean';
        else 
            AgeBin = 'Post-wean';
        end         

        CellCount3 = length(famOverlap3) %used for adding N number to plots 

        figure()
        x = categorical({'FAM vs FAM','FAM vs NOV1','FAM vs NOV2'});
        y = [FAMOverlapC3 DIFFOverlapC3 NOVOverlapC3];%; FAMOverlapC1 NOVOverlapC1; FAMOverlapC4 NOVOverlapC4];
        errors = [ FAMErrC3 DIFFErrC3 NOVErrC3]; %; FAMErrC1  NOVErrC1;FAMErrC4 NOVErrC4];
%          xticks = ({strcat('GCs (',num2str(CellCount2),')'),strcat('MCs (',num2str(CellCount3),')')}); %strcat('INs (',num2str(CellCount1),')'),,strcat('C4(',num2str(CellCount4),')')});
%         errorBars = gra_groupedbars(y, errors);
%         errorBars = set(gca,'xticklabels', xticks);
        RateOver= bar(x,y);
        set(gca, 'FontSize', 16)
        hold on
            er = errorbar(x,y,errors,errors);    
            er.Color = [0 0 0];                            
            er.LineStyle = 'none';  
        hold off
        title(strcat(AgeBin,': P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2)))); %r',num2str(ratBins(itRat,1)),': put this back in when you fix rat binning 
            ylabel('Rate Overlap Score','FontSize', 16)
%             xlabel('Cell Type (N)')
        
%          figure()
%             y = mean(SI_spat(cluster3)) %looks like these were uncommented so I could add means ans erros somewhere
%             errorC3 = transpose(std(SI_spat(cluster3))/sqrt(length(cluster3)));
%             errors = errorC3
%         %         xticks = ({strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2)))}); 
%         %         errorBars = gra_groupedbars(y, errors);
%         %         errorBars = set(gca,'xticklabels', xticks);
%             bar(y)
%             title(strcat('P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2)))); 
%                 ylabel('Spatial Information Score')
%                 xlabel('Age')

%         Age = [Age;itAge*ones(length(famOverlap3)*3,1)];
%         Environment = [Environment;ones(length(famOverlap3),1); 2*ones(length(famOverlap3),1);3*ones(length(famOverlap3),1)];
         
%         %prep parameters for mixed anova 
%         famOverlap3 = atanh(famOverlap3);% to make it parametric for mixed anova
%         novOverlap3 = atanh(novOverlap3);
%         diffOverlap3 = atanh(diffOverlap3);
        
        rateOverlaps = [rateOverlaps;famOverlap3,diffOverlap3,novOverlap3];
        if itAge == 1
            AgeBin1Subjects = length(famOverlap3);
        elseif itAge == 2
            AgeBin2Subjects = length(famOverlap3);
        elseif itAge == 3
            AgeBin3Subjects = length(famOverlap3);       
        end
        %reset values used
        meanRate = spatData.meanRate;
        cluster3 = original_cluster;


    end  


    %ive decided to use simple mixed anova from mathworks 
    datamat = rateOverlaps;
    within_factor_names = {'Environment'};
    between_factors = [ones(AgeBin1Subjects, 1); 2 * ones(AgeBin2Subjects, 1)];% 3 * ones(AgeBin3Subjects, 1)];% ]; %
    between_factor_names = {'Age'};
    [tbl, rm] = simple_mixed_anova(datamat, between_factors, within_factor_names, between_factor_names)


end        

