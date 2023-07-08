function spatialCorrelation_thesis(spatData,clusters)
%SPATIAL CORRELATION remapping outcome measure like
%rate overlap but using spatial overlap. 
%   1.uses map_spatialcorr and rates_transform to make this score 
%   2.organized using age binning and indexes for environment.

    load(spatData,'spatData') 
        meanRate = spatData.meanRate;
        rMap = spatData.rMap;
        rMap1stHalf = spatData.rMap1stHalf;
        rMap2ndHalf = spatData.rMap2ndHalf;
        SI_spat = spatData.SI_spat;
        nSpks = spatData.nSpks;
        
        
    %gather age data from cellInfo 
    
        cellInfo = getCellInfo(spatData);
        age = cellInfo(:,3); 
        rat = cellInfo(:,1); 
        corecorded = cellInfo(:,4);
        
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
                for itTrial = 1: width(spatData.env)
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
    
    
    %load clusters
    load(clusters, 'PCA2_clusters', 'DG_ExCluster','CA3_ExCluster')     
    %make clusters from PCA2_clusters keeping old naming convention for
    %convenience in running the old code. 
    mossy =[];
    granule =[]; %keeping og naming convention for a second to see if this code runs 
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 1
            mossy = [mossy;DG_ExCluster(ii)]; %for WF1_AMR_BI_CCC_new - 1 is mc and 2 is gc 
        elseif PCA2_clusters(ii) == 2
            granule = [granule;DG_ExCluster(ii)];% for WF1_AMR_BI_CCC_r1099 - 1 is gc and 2 is mc 

        end
    end 
    
    cluster =  granule; %CA3_ExCluster; % mossy; %
    clustername =  {'Granule Cells'}; %{'CA3 Pyramidal Cells'}; % {'Mossy Cells'};%
    orange = [1, 0.5, 0];
    color = orange;

    Age =[]; 
    Environment = [];
    spatCorrs = [];

    FAMfieldNumC3 = [];    
    NOVfieldNumC3 = [];    
    DIFFfieldNumC3 =[];
    
    %   age binning 
    
     ageBins   =  [17 20; 21 32]; % list of age bins each spanning from col1:col2
     original_cluster = cluster; 
        
     for itAge=1:size(ageBins,1)
            
            indAge = age>=ageBins(itAge,1) & age<=ageBins(itAge,2); % index for current age bin
    
            %get indices for new clusters - shorten cell cluster to age bin
            %only. 
            ages_indexes_in_spatData = [];
            for it_indAge =1: length(indAge) 
                if true(indAge(it_indAge))
                    ages_indexes_in_spatData = [ages_indexes_in_spatData;it_indAge];
                else
                end
            end
            Newcluster = ismember(cluster,ages_indexes_in_spatData);
            cluster3 = cluster(Newcluster);
            cluster_count = cluster3;
            %remove non-spatial cells from the cluster 
            for itClu = 1: length (cluster3) 
                if any(SI_spat(cluster3(itClu),:) > 0.2) || any(nSpks(cluster3(itClu),1:5) > 100)
                    cluster3(itClu) = cluster3(itClu);
                else
                    cluster3(itClu) = 0;
                end          
            end             
            cluster3 = cluster3(cluster3 ~=0);          

            CellCount3 = length(cluster3) %this is here for a sense check when its running
    
            
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
                fam1map = cell2mat(rMap(cluster3(itRM),FamInd(cluster3(itRM),1)));
                fam2map = cell2mat(rMap(cluster3(itRM),FamInd(cluster3(itRM),2)));
                novmap = cell2mat(NewSquare(cluster3(itRM),1));
                diffmap = cell2mat(rMap(cluster3(itRM),DiffInd(cluster3(itRM),1)));

                FAMSpatCorrC3(itRM, :) = manual_map_spatialcorr(fam1map, fam2map);
                NOVSpatCorrC3(itRM, :) = manual_map_spatialcorr(fam2map, novmap);
                SELFSpatCorrC3(itRM, :) = manual_map_spatialcorr(cell2mat(rMap1stHalf(cluster3(itRM),FamInd(cluster3(itRM),1))), cell2mat(rMap2ndHalf(cluster3(itRM),FamInd(cluster3(itRM),1))));
                DIFFSpatCorrC3(itRM,:) = manual_map_spatialcorr(fam2map,diffmap);

            end 
    
            %make means  
            cells_included = length(FAMSpatCorrC3);
            FAMOverlapC3 = mean(rmmissing(FAMSpatCorrC3(FAMSpatCorrC3~=0)));
            NOVOverlapC3 = mean(rmmissing(NOVSpatCorrC3(NOVSpatCorrC3~=0)));    
            DIFFOverlapC3 = mean(rmmissing(DIFFSpatCorrC3(DIFFSpatCorrC3~=0))); 
            SELFOverlapC3 = mean(rmmissing(SELFSpatCorrC3(SELFSpatCorrC3~=0)));  
            % make errors
            FAMErrC3 = std(rmmissing(FAMSpatCorrC3))/sqrt(length(cluster3)); 
            NOVErrC3 = std(rmmissing(NOVSpatCorrC3))/sqrt(length(cluster3));    
            DIFFErrC3 = std(rmmissing(DIFFSpatCorrC3))/sqrt(length(cluster3));
            SELFErrC3 = std(rmmissing(SELFSpatCorrC3))/sqrt(length(cluster3));
    
            %make figures

            %make Agebin titles 

            if itAge == 1
                AgeBin = 'Pre-wean';
            else 
                AgeBin = 'Post-wean';
            end

            %plot of spatial correlations for current age bin/cell type
            figure ()
            x = categorical({'FAM vs FAM','FAM vs NOV1','FAM vs NOV2'});
            y = [FAMOverlapC3 DIFFOverlapC3 NOVOverlapC3];%;FAMOverlapC1 NOVOverlapC1; FAMOverlapC4 NOVOverlapC4];
            xticks = ({strcat('MCs : ',num2str(CellCount3))});%,strcat('INs : ',num2str(CellCount1)),strcat('C4: ',num2str(CellCount4))});
            errors = [FAMErrC3 DIFFErrC3 NOVErrC3];%; FAMErrC1  NOVErrC1;FAMErrC4 NOVErrC4];
    %         errorBars = gra_groupedbars(y, errors);
            SpatialCorr= bar(x,y);
            set(gca, 'FontSize', 16)
            hold on
                er = errorbar(x,y,errors,errors);    
                er.Color = [0 0 0];                            
                er.LineStyle = 'none';  
            hold off
    %         errorBars = set(gca,'xticklabels', xticks);
            title(strcat(AgeBin,': P',num2str(ageBins(itAge,1)),' to P',num2str(ageBins(itAge,2))));
                    ylabel('Spatial Correlation', 'FontSize', 16);
%                     xlabel('Cell Type: N');
                    ylim([0 0.6]);

            
            %plot of number of fields for current age bin/cell type
            [field_num] = field_counter(rMap, meanRate, SELFSpatCorrC3, cluster3);

%              figure; 
%              histogram(field_num,'FaceColor', color, 'FaceAlpha', 0.5);
%              set(gca, 'FontSize', 16) % Adjust the font size to your preference
% %              titlename = ['Number of Fields for' AgeBin ': ' clustername];
%              titlename = ['Number of Fields for: ' clustername];
%              title(titlename, 'FontSize', 14) % Adjust the font size to your preference
%              xlabel('Number of Fields', 'FontSize', 16) % Adjust the font size to your preference
%              ylabel('Cell Number', 'FontSize', 16) % Adjust the font size to your preference
%              percentage_of_cells_with_fields = length(field_num)/length(cluster_count)*100

             % stats
             FAMSpatCorrC3 = atanh(FAMSpatCorrC3);% to make it parametric for mixed anova
             NOVSpatCorrC3 = atanh(NOVSpatCorrC3);
             DIFFSpatCorrC3 = atanh(DIFFSpatCorrC3);
    
             %factors for n-way anova 
    %          Environment = [Environment;ones(length( NOVSpatCorrC3),1); 2*ones(length( NOVSpatCorrC3),1);3*ones(length( NOVSpatCorrC3),1)];        
    %          Age = [Age;itAge*ones(length(FAMSpatCorrC3)*3,1)];%;repmat(ageBins(itAge,1),length(cluster3),1);repmat(ageBins(itAge,1),length(cluster3),1)];% *2 removed for KW atempt 
    %          spatCorrs = [spatCorrs;FAMSpatCorrC3 ;DIFFSpatCorrC3; NOVSpatCorrC3];%[rateOverlaps; novOverlap2 - famOverlap2];%[famOverlap2 ;novOverlap2];%[rateOverlaps;famOverlap3 ,novOverlap3]%
             
             %factors for mixed anova 
             spatCorrs = [spatCorrs;FAMSpatCorrC3 ,DIFFSpatCorrC3, NOVSpatCorrC3];
             if itAge == 1
                 AgeBin1Subjects = length(FAMSpatCorrC3);
             elseif itAge == 2
                AgeBin2Subjects = length(FAMSpatCorrC3);
             elseif itAge == 3
                AgeBin3Subjects = length(FAMSpatCorrC3);       
             end
    
             %reset variables
             CellCount3 = 0;
             cluster3= original_cluster;
            
     end  

    % mixed anova 

    datamat = spatCorrs;
    within_factor_names = {'Environment'};
    between_factors = [ones(AgeBin1Subjects, 1); 2 * ones(AgeBin2Subjects, 1)];%; 3 * ones(AgeBin3Subjects, 1)];
    between_factor_names = {'Age'};
    [tbl, rm] = simple_mixed_anova(datamat, between_factors, within_factor_names, between_factor_names)
   
    % n - way anova 
    %     %t-test for spatial score and age
    %        %[p,tbl] = anova1(SI_score,Ages);
    %        [p,tbl,stats]=anovan(spatCorrs,{Age,Environment},'model','interaction', 'varnames',{'Age','Environment'});
    %         c = multcompare(stats)
    % %     save('spatRemData.mat','gc_spatRem','mc_spatRem')

end 

function [r] = manual_map_spatialcorr(rates1, rates2, varargin)

    if isempty(rates1) || isempty(rates2)
        r = NaN;  return
    end
    if max(rates1(:))==0 || max(rates2(:))==0
        r = NaN;  return
    end
    
    rates1 = double(rates1);   
    rates2 = double(rates2);
    
    visInd = ~isnan(rates1) & ~isnan(rates2);   % Find unvisited bins in either trial.
    
    if ~isempty(varargin) 
        if strcmp(varargin{1},'ignoreMutualZero')
            zeroInd = rates1==0 & rates2==0;
            visInd = visInd & ~zeroInd;
        else
            error('Input argument not recognized');
        end
    end
    
    % Manual calculation
    r_mat = corrcoef(rates1(visInd), rates2(visInd));
    r = r_mat(1, 2);
end


function [field_num] = field_counter(rMap, meanRate,SELFSpatCorrC3, curr_cluster)
%make number of fields to be used in plots 
        field_num = [];
        fieldCount = 0;
        for itFC = 1: length (curr_cluster)            
            [~, max_rate_ind] = nanmax(meanRate(curr_cluster(itFC),1:5));
            map_rate_thr = cell2mat(rMap(curr_cluster(itFC),max_rate_ind)) >= meanRate(curr_cluster(itFC),max_rate_ind);
            stats = regionprops(map_rate_thr, 'Area');
            Areas = cell2mat(struct2cell(stats));
            for itA = 1:length(Areas)
                if Areas(itA) >= 25 && Areas(itA) <= 150
                    fieldCount = fieldCount + 1;
                end
            end
             if fieldCount ~= 0 %&& SELFSpatCorrC3(itFC) > 0.05
                field_num =[field_num; fieldCount];
             else
                 field_num = field_num;
             end
            fieldCount = 0;
        end
end
