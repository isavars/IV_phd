function activeCells()

    load ('r889_P18_spatData.mat', 'spatData')
    
    nSpks = spatData.nSpks; %load spike numbers from file 
    env = spatData.env; %load environments 
    corecordedNo = corecordedCells();
    
    nSpks(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    
    corecorded = zeros(1, 5);
    
%     for jj = 1:length(corecordedNo)
%         corecorded (jj,:) = corecordedNo(jj)
%     end

    for jj = 1:length(nSpks)
        corecorded(jj,:)= corecordedNo(jj);
        if nSpks(jj,:) > 10
            corecorded(jj)= corecordedNo(jj);
        elseif isnan (nSpks (jj,:))
            corecordedNo (jj) = NaN;
        else nSpks(jj,:) <= 10
            corecorded(jj)= corecordedNo(jj) -1;
        end
    end
    corecordedNo
    whos
end