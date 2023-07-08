function cell_property_proportions (data,clusters)
%this is going to make numbers that go into pie charts made in illustrator
%- the pie charts represent proportions of cell populations that do a
%certain thing. This could just be a subfunction in rateOverlap_thesis
%if it uses all the same parameters and has the same loop logic (by cluster
%and by age...) 

    %load spatial data for all cells and gather necessary parameters 
    load(data, 'spatData')
    meanRate = spatData.meanRate;

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


    %load and make clustersfor each cell type 
    load(clusters, 'PCA2_clusters', 'DG_ExCluster','CA3_ExCluster')     
    %make clusters from PCA2_clusters keeping old naming convention for
    %convenience in running the old code. 
    mossy_cluster =[];
    granule_cluster =[]; %keeping og naming convention for a second to see if this code runs 
    for ii = 1: length(PCA2_clusters)
        if PCA2_clusters(ii) == 1
            granule_cluster = [granule_cluster;DG_ExCluster(ii)]; %for WF1_AMR_BI_CCC_new - 2 is gc and 1 is mc 
        elseif PCA2_clusters(ii) == 2
            mossy_cluster = [mossy_cluster;DG_ExCluster(ii)];% for WF1_AMR_BI_DS2_new - 2 is gc and 1 is mc 
        end
    end
    pyramidal_cluster = CA3_ExCluster;

    clusters = {granule_cluster, mossy_cluster, pyramidal_cluster}; 
    clustername = {'granule', 'mossy', 'CA3'}; %to me used as tile change also depending on clusters plotted. 

    %loop over spatData and assign excitatory cells to age bins

    %gather age data from cellInfo 
    cellInfo = getCellInfo(spatData);
    age = cellInfo(:,3);

    %separate the data into age bins 
    
    ageBins = [17 20; 21 31]; %21 24; 
 
        
    for ii = 1: length(clusters) 
        cluster = clusters{ii};
        original_cluster = cluster;
        for itAge = 1: size(ageBins,1)            
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
            cluster = cluster(Newcluster);
            %bring in necessary values 
            meanRate_clu = spatData.meanRate(cluster,:);
            nSpks_clu = spatData.nSpks(cluster,:);
            SI_spat_clu = spatData.SI_spat(cluster,:);
            FamInd_clu = FamInd(cluster,:);
            NovInd_clu = NovInd(cluster);
            DiffInd_clu = DiffInd(cluster); 

            %make environment indexes 

            %make proportion of active cells per envirnoment
            active_in_sleep_only = 0;
            active_in_1_environment = 0;
            active_in_2_environments = 0;
            active_in_all_environments = 0; 
            %proportion of cells that are spatial (pick critera - should be same for using spatcorr)
            is_spatial = 0;
            is_non_spatial = 0; 
            for jj = 1: length(cluster)
                %decide wht makes a cell active? also this needs to deal with
                %environment names to deal with fam box being 3 trials - change
                %in the mornig 
                if all(nSpks_clu(jj,1:5) < 100)
                    active_in_sleep_only = active_in_sleep_only +1;
                elseif (sum(nSpks_clu(jj,FamInd_clu(jj,1):FamInd_clu(jj,2)) > 100) == 2) + sum(nSpks_clu(jj,NovInd_clu(jj)) > 100)  + sum(nSpks_clu(jj,DiffInd_clu(jj)) > 100) == 1
                    active_in_1_environment = active_in_1_environment +1;
                elseif (sum(nSpks_clu(jj,FamInd_clu(jj,1):FamInd_clu(jj,2)) > 100) == 2) + sum(nSpks_clu(jj,NovInd_clu(jj)) > 100)  + sum(nSpks_clu(jj,DiffInd_clu(jj)) > 100) == 2
                    active_in_2_environments = active_in_2_environments +1;  
                elseif (sum(nSpks_clu(jj,FamInd_clu(jj,1):FamInd_clu(jj,2)) > 100) == 2) + sum(nSpks_clu(jj,NovInd_clu(jj)) > 100)  + sum(nSpks_clu(jj,DiffInd_clu(jj)) > 100) == 3
                    active_in_all_environments = active_in_all_environments +1; 
                else 
                    active_in_sleep_only = active_in_sleep_only +1;
                end
                
                if any(SI_spat_clu(jj,:) > 0.5) && any(nSpks_clu(jj,1:5) > 100)
                    is_spatial = is_spatial +1;
                else
                    is_non_spatial = is_non_spatial +1;
                end 
            end        
%             display(['Percentages for ' clustername{ii} ':' ]);
            active_in_sleep_only = active_in_sleep_only./length(cluster)*100;
            active_in_1_environment= active_in_1_environment./length(cluster)*100;
            active_in_2_environments = active_in_2_environments./length(cluster)*100;
            active_in_all_environments = active_in_all_environments./length(cluster)*100;
            
            is_spatial = is_spatial./length(cluster)*100;
            is_non_spatial = is_non_spatial./length(cluster)*100;

            %make Agebin titles 

            if itAge == 1
                AgeBin = 'Pre-wean';
            else 
                AgeBin = 'Post-wean';
            end
                
            %make pie charts
            values = [active_in_sleep_only, active_in_1_environment, active_in_2_environments, active_in_all_environments];
            lables = {'active in sleep only', 'active in 1 environment', 'active in 2 environments', 'active in all environments'};
%             figure;
%             pie(values,lables);
%             title(['Proportion of cells active in each environment for: ' clustername{ii} ' (' AgeBin ')'])

            values = [is_spatial, is_non_spatial];
            lables = {'is spatial', 'is non-spatial'};
            figure;
            pie(values,lables);
            title(['Proportion of cells spatial vs non-spatial for: ' clustername{ii} ' (' AgeBin ')'])

            %reset variables 
            cluster =original_cluster;
        end
        
        
        

    end 
end 