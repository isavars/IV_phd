function spatialCorrelation_tetrodes_thesis(spatData,clusters)
%SPATIAL CORRELATION remapping outcome measure like
%rate overlap but using spatial overlap. 
%   1.uses map_spatialcorr and rates_transform to make this score 
%   2.organized using age binning and indexes for environment.

%needs to be adaptable for tetrodes 

    load(spatData,'spatData') 
        meanRate = spatData.meanRate;
        peakRate = spatData.peakRate;
        rMap = spatData.rMap;
        rMap1stHalf = spatData.rMap1stHalf;
        rMap2ndHalf = spatData.rMap2ndHalf;
        SI_spat = spatData.SI_spat;
        nSpks = spatData.nSpks;
                
    %gather age data from cellInfo 
    
        cellInfo = getCellInfo(spatData);
        age = cellInfo(:,3); 
    
    %make indexes for sleep and wake trials 
    sleepMeanRate_all= zeros(size(spatData,1),1);
    sleep_idx = zeros(size(spatData,1),1);
    awakeMeanRate_all = zeros(size(spatData,1),1);
    maxAwakeMeanRate_all = zeros(size(spatData,1),1);
    wake_idx = cell(size(spatData,1),1);
    for itCl = 1: height(spatData)
        sleep_trials = strcmp(string(spatData.env(itCl,:)),'sleep');
        sleep_idx_temp = find(sleep_trials);
        if size(sleep_idx_temp,2) > 1 
            sleep_idx(itCl) = sleep_idx_temp(2); %dealing with trial with more than one sleep
        else 
            sleep_idx(itCl) = sleep_idx_temp;
        end
        nov_trials = strcmp(string(spatData.env(itCl,:)),'nov');
        fam_trials = strcmp(string(spatData.env(itCl,:)),'fam');
        wake_trials = nov_trials + fam_trials; 
        %datasets have different numbers of wake trials 
        wake_idx_temp = find(wake_trials);
        wake_idx{itCl} = wake_idx_temp;
        awakeMeanRate_all(itCl) = nanmean(spatData.meanRate(itCl,wake_idx_temp));
        maxAwakeMeanRate_all(itCl) = nanmax(spatData.meanRate(itCl,wake_idx_temp));
        sleepMeanRate_all(itCl) = nanmean(spatData.peakRate(itCl,sleep_idx_temp));
    end 
        
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
                        %add top firing rate trials 
                        FamIndT = [FamIndT, itTrial];
                        FamIndT = nonzeros(FamIndT)';
                        famCount = famCount + 1;
                        if famCount > maxFam                            
                            nSpks_per_fam_trial = nSpks(itCell,FamIndT);
                            [~, idx] = sort(nSpks_per_fam_trial, 'descend'); %find top 2 firing rate trials 
                            top_two = FamIndT(idx(1:2));
                            FamInd(itCell,:)= top_two;
                        else  
                            FamInd(itCell,famCount)= FamIndT(famCount);
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
                FamIndT = [];
            end
    
    
    %load clusters
    load(clusters, 'PCA2_clusters', 'DG_ExCluster','CA3_ExCluster')     
    %make clusters from PCA2_clusters keeping old naming convention for
    %convenience in running the old code. 
    mossy =[];
    granule =[]; %keeping og naming convention for a second to see if this code runs 
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 2
            mossy = [mossy;DG_ExCluster(ii)]; 
        elseif PCA2_clusters(ii) == 1
            granule = [granule;DG_ExCluster(ii)];
        end
    end 
    
    clusters =  {granule, mossy}; %granule; %CA3_ExCluster; % 
    clusternames =  {'Granule Cells','Mossy Cells'};%{'CA3 Pyramidal Cells'}; % 
    orange = [1, 0.5, 0];
    colors = {orange, 'g'};

    for itC = 1: length(clusters)

    cluster = clusters{itC};
    clustername = clusternames{itC};
    color1 = colors{itC};

        Age =[]; 
        Environment = [];
        spatCorrs = [];
    
        FAMfieldNumC3 = [];    
        NOVfieldNumC3 = [];    
        DIFFfieldNumC3 =[];
        
        %   age binning 
        
         ageBins   =  [17 20; 21 32; 40 40]; % list of age bins each spanning from col1:col2
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
                    if any(spatData.sig_SI(cluster3(itClu),wake_idx{cluster3(itClu)}) == 1) && any(nSpks(cluster3(itClu),wake_idx{cluster3(itClu)}) > 100) && any(SI_spat(cluster3(itClu),wake_idx{cluster3(itClu)}) > 0.2) %&& any(nSpks(cluster3(itClu),wake_idx{cluster3(itClu)}) > 75)%&& any(SI_spat(cluster3(itClu),1:5) > 0.2) %|| any(nSpks(cluster3(itClu),1:5) > 100) %any(SI_spat(cluster3(itClu),1:5) > 0.5) && any(nSpks(cluster3(itClu),1:5) > 75) && any(spatData.sig(cluster3(itClu)) == 1) %these are the knierm filters 
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
    
                %deal with tetrode trials 
                if isempty(DiffInd)
                    %display('no diff trial, ignore middle column')
                    DiffInd = FamInd;
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
                x = categorical({'FAM vs FAM','FAM vs NOV2'}); %'FAM vs NOV1',
                y = [FAMOverlapC3  NOVOverlapC3];%;FAMOverlapC1 NOVOverlapC1; FAMOverlapC4 NOVOverlapC4]; DIFFOverlapC3
                xticks = ({strcat('MCs : ',num2str(CellCount3))});%,strcat('INs : ',num2str(CellCount1)),strcat('C4: ',num2str(CellCount4))});
                errors = [FAMErrC3 NOVErrC3];%; FAMErrC1  NOVErrC1;FAMErrC4 NOVErrC4]; DIFFErrC3
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

                %calculat peak rates and SI scores of spatial cells for
                %stats reporting 

                SI_scores = [];
                peak_rates = [];
                for ii = 1: length(cluster3)
                    [~, maxSpksPos] = nanmax(spatData.nSpks(cluster3(ii),wake_idx{cluster3(ii)})); 
                    SI_spat_temp = spatData.SI_spat(cluster3(ii),maxSpksPos);
                    SI_scores = [SI_scores; SI_spat_temp];
                    peak_rates_temp = spatData.peakRate(cluster3(ii),maxSpksPos);                    
                    peak_rates = [peak_rates; peak_rates_temp];
                end    
                
                %calculate and plot of number of fields for current age bin/cell type
                [field_num, singlefield_count,multifield_count] = field_counter(rMap, peakRate, meanRate,cluster3, color1, clustername, AgeBin, cluster_count, wake_idx);

                %calculate proporiton of cells with fields (spatial vs non
                %spatial
                cells_with_fields = sum(field_num ~= 0);
                cells_without_fields = length(cluster_count)-sum(field_num ~= 0);

                %make pies for spatial vs non-spatial
                spatiality_values= [cells_with_fields,cells_without_fields] ;
                spatial_cells = length(field_num);    
                nonSpatial_cells = length(cluster_count) - spatial_cells;
                field_count_all_cells = [field_num; zeros(nonSpatial_cells,1)];                
%                 spatiality_values=[spatial_cells,nonSpatial_cells ];
                spatiality_lables = {'Spatial ', 'Non-spatial '};
%                 make_pies(spatiality_values,spatiality_lables, clustername, AgeBin);



                %calculate proportion of multifield vs single field 
                percentage_of_multifield_cells = (multifield_count/ (multifield_count + singlefield_count))*100;
                percentage_of_singlefield_cells = (singlefield_count/ (multifield_count + singlefield_count))*100;
                %make pies for feild number
                fieldness_values= [multifield_count,singlefield_count] ;
                fieldness_lables = {'Multi field ', 'Single field'};
                make_pies(fieldness_values,fieldness_lables, clustername, AgeBin);
  
                %calculate silent vs active proportions (this one can also
                %be for multiple envs remapping results)
                
                [percentage_of_active_cells, percentage_of_silent_cells, percentage_active_in_all, percentage_active_in_one, percentage_active_in_fam_and_diff, percentage_active_in_fam_and_nov] = silent_vs_active(nSpks, wake_idx, sleep_idx, cluster_count, FamInd, DiffInd, NovInd, maxAwakeMeanRate_all, sleepMeanRate_all, spatData);
                 %make pies for silent vs active 
                activity_values= [percentage_of_silent_cells,percentage_of_active_cells] ;
                activity_lables = {'Silent ', 'Active '};
                make_pies(activity_values,activity_lables, clustername, AgeBin);

                remapping_values = [percentage_active_in_all, percentage_active_in_one, percentage_active_in_fam_and_diff, percentage_active_in_fam_and_nov];
                remapping_labels = {'Active in all', 'Active in one Base', 'Active in Base only', 'Active in Diff only'};
                % need to prep data for stats on proporitons
                make_pies(remapping_values,remapping_labels, clustername, AgeBin);
                
                %get values for proporitons statisitics %chaneg values
                %depending on test you want to run 

                values = fieldness_values;
                value_name = 'fieldness_values';

                if itC == 1
                    if itAge == 1
                        granule_prewean = values;
                    elseif itAge ==2
                        granule_postwean = values;
                    elseif itAge ==3
                        granule_adult = values;
                    end 
                elseif itC == 2
                    if itAge == 1
                        mossy_prewean = values;
                    elseif itAge ==2
                        mossy_postwean = values;
                    elseif itAge ==3
                        mossy_adult = values;                        
                    end 
                end                
               
                 % stats
                 FAMSpatCorrC3 = atanh(FAMSpatCorrC3);% to make it parametric for mixed anova
                 NOVSpatCorrC3 = atanh(NOVSpatCorrC3);
                 %DIFFSpatCorrC3 = atanh(DIFFSpatCorrC3);
        
                 %factors for mixed anova 
                 %spatial correlations
                 spatCorrs = [spatCorrs;FAMSpatCorrC3 , NOVSpatCorrC3];%[spatCorrs;FAMSpatCorrC3 ,DIFFSpatCorrC3, NOVSpatCorrC3];
                 %age bins
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

        %% STATS
    
        % mixed anova 
    
        datamat = spatCorrs;
        within_factor_names = {'Environment'};
        between_factors = [ones(AgeBin1Subjects, 1); 2 * ones(AgeBin2Subjects, 1); 3 * ones(AgeBin3Subjects, 1)];
        between_factor_names = {'Age'};
        [tbl, rm] = simple_mixed_anova(datamat, between_factors, within_factor_names, between_factor_names)
        
        %save table for spss - so you cna do the post hoc 
        % Assuming spatCorrs is of size [numSubjects x numEnvironments]
        [numSubjects, numEnvironments] = size(spatCorrs);
        
        % Create a table for export
        SubjectID = repmat((1:numSubjects)', numEnvironments, 1);
        Environment = repelem((1:numEnvironments)', numSubjects, 1);
        
        % Convert the matrix to long format
        SpatialCorrelation = reshape(spatCorrs, [], 1);
        
        % If AgeBin1Subjects, AgeBin2Subjects are the counts of subjects in each age bin
        Ages = [ones(AgeBin1Subjects, 1); 2* ones(AgeBin2Subjects, 1); 3 * ones(AgeBin3Subjects, 1)];
        Age = [Ages;Ages]; %there's probably a better way to do that but who cares - repreating ages per num of environments
        
        % Combine everything into a table
        T = table(SubjectID, Age, Environment, SpatialCorrelation);
        writetable(T, [clustername '_tetrodes_spatCorrs_for_SPSS.csv'])

    end

    %table for SPSS input - to compare proportions - remove CA3 cluster
    % Create an array to store the data
    data = [granule_prewean; mossy_prewean; granule_postwean; mossy_postwean; granule_adult; mossy_adult];

    % Create Age and CellType arrays
    Age = [repmat({'prewean'}, [2, 1]); repmat({'postwean'}, [2, 1]); repmat({'adult'}, [2, 1])];
    CellType = {'granule'; 'mossy'; 'granule'; 'mossy'; 'granule'; 'mossy'};
    
    % Convert data to table
    T = table(Age, CellType, data(:,1), data(:,2), 'VariableNames', {'Age', 'CellType', 'Spatial', 'NonSpatial'});
    
    % Display the table
    disp(T);

%     %table for spss when comparing individual values per cell in an anova 
% 
%     % Create data array
%     data = [granule_prewean; mossy_prewean; granule_postwean; mossy_postwean; granule_adult; mossy_adult];
%     
%     % Create Age and CellType arrays
%     Age = [repmat({'prewean'}, [length(granule_prewean) + length(mossy_prewean),1 ]); 
%            repmat({'postwean'}, [length(granule_postwean) + length(mossy_postwean), 1]);
%            repmat({'adult'}, [length(granule_adult) + length(mossy_adult), 1])];
%            
%     CellType = [repmat({'granule'}, [length(granule_prewean), 1]);
%                 repmat({'mossy'}, [length(mossy_prewean), 1]);
%                 repmat({'granule'}, [length(granule_postwean), 1]);
%                 repmat({'mossy'}, [length(mossy_postwean), 1]);
%                 repmat({'granule'}, [length(granule_adult), 1]);
%                 repmat({'mossy'}, [length(mossy_adult), 1])];
% 
%     disp(length(Age));
%     disp(length(CellType));
%     disp(length(data));
% 
%     
%     % Convert data to table
%     T = table(Age, CellType, data, 'VariableNames', {'Age', 'CellType', value_name});
%     
%     % Display the table
%     disp(T);
% %     swtest(data)
% 
%     writetable(T, [value_name '_tetrode_data_for_spss.csv']);
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


function [field_num, singlefield_count,multifield_count] = field_counter(rMap, peakRate, meanRate, curr_cluster, color1, clustername, AgeBin, cluster_count, wake_idx)
%make number of fields to be used in plots 
%needs to grab the trial with the most fields/real fields -check knierim
%-using highest firing rate 
%add cells with zero feilds to the count 
        
        num_non_spatial_cells = length(cluster_count) - length(curr_cluster);
        
        field_num = [];
        fieldCount = 0;
        for itFC = 1: length (curr_cluster)            
            [~, max_rate_ind] = nanmax(peakRate(curr_cluster(itFC),wake_idx{curr_cluster(itFC)})); %using the trial with the highest mean rate
            map_rate_thr = cell2mat(rMap(curr_cluster(itFC),max_rate_ind)) >= peakRate(curr_cluster(itFC),max_rate_ind)*0.3;
            stats = regionprops(map_rate_thr, 'Area');
            Areas = cell2mat(struct2cell(stats));
            for itA = 1:length(Areas)
                if Areas(itA) >= 25 && Areas(itA) <= 200 %changed min to 15
                    fieldCount = fieldCount + 1;
                end
            end

            field_num =[field_num; fieldCount];

            fieldCount = 0;
        end

%         field_num = [field_num; zeros(num_non_spatial_cells,1)];
        
        %make multi and single field counts 
        multifield_count = 0;
        singlefield_count = 0;
        for ii = 1: length(field_num)
            if field_num(ii) == 1
                singlefield_count = singlefield_count + 1;
            elseif field_num(ii) > 1
                multifield_count = multifield_count + 1;
            end 
        end 
        
        spatial_cells = length(field_num);    
        nonSpatial_cells = length(cluster_count) - spatial_cells;
        field_count_all_cells = [field_num; zeros(nonSpatial_cells,1)]; 
        %plot field numbers 
         figure; 
         histogram(field_count_all_cells,'FaceColor', color1, 'FaceAlpha', 0.5);
         set(gca, 'FontSize', 16) % Adjust the font size to your preference
         titlename = ['Number of Fields for ' AgeBin ': ' clustername];
%          titlename = ['Number of Fields for: ' clustername];
         title(titlename, 'FontSize', 14) % Adjust the font size to your preference
         xlabel('Number of Fields', 'FontSize', 16) % Adjust the font size to your preference
         ylabel('Cell Number', 'FontSize', 16) % Adjust the font size to your preference 


end

function make_pies(values,lables, clustername, AgeBin)
    %takes in arrays of data and produces pie charts with the correct
    %labels maybe this should also be called for the stats 
    figure;
    pie(values,lables);
    title(['Proportion of ' lables{1} ' vs ' lables{2} ' for: ' clustername ' (' AgeBin ')'])
end

function [active_count, silent_count, percentage_active_in_all, percentage_active_in_one, percentage_active_in_fam_and_diff, percentage_active_in_fam_and_nov] = silent_vs_active(nSpks, wake_idx, sleep_idx, cluster_count, FamInd, DiffInd, NovInd, maxAwakeMeanRate_all, sleepMeanRate_all, spatData)
    %make active cells and silent cell proporitons - could also do active
    %in different environments for this one for the remapping chapter -
    %after discussing with tom what this should look like
    active_count= 0;
    silent_count = 0;
    active_in_all = 0;
    active_in_one = 0;
    active_in_fam_and_diff = 0;
    active_in_fam_and_nov = 0;

    for ii = cluster_count'
        if ~any(nSpks(ii,wake_idx{ii}) > 100) && (maxAwakeMeanRate_all(ii) < sleepMeanRate_all(ii))%changed to 100 make it a bit stricter
            silent_count = silent_count +1;
        else
            active_count = active_count +1;
        end
%         if any(nSpks(ii,wake_idx{ii}) > 100) %&& (awakeMeanRate_all(ii) > sleepMeanRate_all(ii))%changed to 100 make it a bit stricter
%             active_count = active_count +1;
%         else
%             silent_count = silent_count +1;
%         end

        if all(nSpks(ii,wake_idx{ii}) > 100)
            active_in_all = active_in_all + 1;
        elseif nSpks(ii,FamInd(ii,1)) > 100 && nSpks(ii,FamInd(ii,2)) > 100
            active_in_fam_and_diff = active_in_fam_and_diff + 1; %changing this one to active in all baselines name stays the same but postion and label are good
        elseif nSpks(ii,NovInd(ii)) > 100
            active_in_fam_and_nov = active_in_fam_and_nov + 1; %and this one is active in nov only 
        elseif sum([nSpks(ii,FamInd(ii,1)) > 100, nSpks(ii,FamInd(ii,2)) > 100, nSpks(ii,NovInd(ii)) > 100]) == 1
            active_in_one = active_in_one + 1;%this is a random baseline but no other tirials 
        end


    end

%     active_count = 0;
%     silent_count = 0;
%     active_in_all = 0;
%     active_in_one = 0;
%     active_in_fam_and_diff = 0;
%     active_in_fam_and_nov = 0;
%     
%     isActiveOverall = false(length(cluster_count),1);  
%     isActivePerEnv = false(length(cluster_count), 5);  % Initialize a matrix to check activity per environment.
% 
% 
%     clust_idx = 0;
%     
%     % First, determine overall activity of the cell.
%     for ii = cluster_count'
%         clust_idx = clust_idx +1;
%         % Identifying silent vs active cell
%         if ~any(nSpks(ii,wake_idx{ii}) > 100) && (maxAwakeMeanRate_all(ii) < sleepMeanRate_all(ii))
%             silent_count = silent_count + 1;
%         else
%             active_count = active_count + 1;
%             isActiveOverall(clust_idx,:) = true;
%         end
%     end
%     
%     clust_idx = 0;
%     % If a cell is active overall, check its activity per environment.
%     for ii = cluster_count'
%         clust_idx = clust_idx +1;
%         if isActiveOverall(clust_idx)
%             for env = 1:length(wake_idx{ii})
%                 if nSpks(ii, wake_idx{ii}(env)) > 100 && spatData.meanRate(ii, wake_idx{ii}(env)) >= sleepMeanRate_all(ii)
%                     isActivePerEnv(clust_idx, env) = true;
%                 end
%             end
%     
%             if all(isActivePerEnv(clust_idx, :))
%                 active_in_all = active_in_all + 1;
%             end
%             if sum([isActivePerEnv(clust_idx, FamInd(ii)), isActivePerEnv(clust_idx, DiffInd(ii)), isActivePerEnv(clust_idx, NovInd(ii))]) == 1
%                 active_in_one = active_in_one + 1;
%             end
%             if isActivePerEnv(clust_idx, FamInd(ii,1)) && isActivePerEnv(clust_idx, DiffInd(ii))
%                 active_in_fam_and_diff = active_in_fam_and_diff + 1;
%             end
%             if isActivePerEnv(clust_idx, FamInd(ii,1)) && isActivePerEnv(clust_idx, NovInd(ii))
%                 active_in_fam_and_nov = active_in_fam_and_nov + 1;
%             end
%         end
%     end


    percentage_of_active_cells = active_count/length(cluster_count) *100;
    percentage_of_silent_cells = silent_count/length(cluster_count) *100;

    percentage_active_in_all = active_in_all/length(cluster_count) *100;
    percentage_active_in_one = active_in_one/length(cluster_count) *100;
    percentage_active_in_fam_and_diff = active_in_fam_and_diff/length(cluster_count) *100;
    percentage_active_in_fam_and_nov = active_in_fam_and_nov/length(cluster_count) *100;


end 

