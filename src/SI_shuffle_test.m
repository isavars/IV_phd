function spatData = SI_shuffle_test(data, shuffles)
%itterate over each line of spatData and compare SI scores to 100 SI scores
%produced by shuffling spike times and positons. Gather SI_shuffles and for
%each cell thats on both tables create a new varible entry in spatData
%sig_SI - this code should work but try on the final shuffles. 

%parameters 
maxTrials = 6;

load(data,'spatData')    
load(shuffles, 'SI_shuffles')

%data collection + appending 
    for itC = 1: height(spatData) 
        for itT = 1: maxTrials
            valueToFind = spatData.cellID{itC};
            matching_row = [];
            for i = 1:length(SI_shuffles.cellID)
                if strcmp(SI_shuffles.cellID{i}, valueToFind)
                    matching_row = i;
                end
            end
            if ~isempty(matching_row)
                originalScore = spatData.SI_spat(itC,itT);
                shuffledScores = SI_shuffles.SI{matching_row,itT};
                percentage = 90; % for 95%, use 90 for 90%, etc.

                thresholdScore = prctile(shuffledScores, percentage);

                %test if shuffled scores are less than the original score
                if originalScore > thresholdScore %all(shuffledScores < originalScore) %do different percentiles of the og score (95th,90th and 75th)
                    spatData.sig_SI(itC,itT) = true;
                else
                    spatData.sig_SI(itC,itT) = false;
                end
            end
        end
    end 
end