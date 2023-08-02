function SI_shuffle_test(data, shuffles)
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
        count = 0;
        for itT = 1: maxTrials
            valueToFind = spatData.cellID{itC};
            matching_row = [];
            for i = 1:length(SI_shuffles.cellID)
                if strcmp(SI_shuffles.cellID{i}, valueToFind)
                    matching_row = [matching_row, i];
                end
            end
            if ~isempty(matching_row)
                originalScore = spatData.SI_spat{itC,itT};
                shuffledScores = SI_shuffles.SI{matching_row,itT};
                count = count +1;
                %test if all shuffled scores are less than the original score
                if all(shuffledScores < originalScore)
                    spatData.sig_SI(itC) = true;
                else
                    spatData.sig_SI(itC) = false;
                end
            end
        end
    end 
end