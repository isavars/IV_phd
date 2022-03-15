function corecorded = corecordedCellsV2(spatData)
%find the number of co-recorded cells per tetrode 
% use tetNo and date to count co-recorded cells per tetrode 
% Also- make a cell # by environment array that shows how many cells were
% corecorded in each environment 
    load ('r889_P18-P22_spatData.mat', 'spatData')
    
    cellID = spatData.cellID;
    dataset = spatData.dataset;
%     nSpks = spatData.nSpks; %load spike numbers from file 
%     env = spatData.env; %load environments 
%     nSpks(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    
    tetNo_ID= zeros(length(cellID), 2);%array to contain cell ID info including co-recorded cells 
    
    for jj= 1: length(cellID)
        newCellID= string(cellID{jj});
        tetNo_ID(jj,2)= double(extractBetween(newCellID,'t','c')); %tetrode number
    end
    
    output = zeros(length(cellID), 3);
    counter=0;
    unique_IDs = unique(dataset)
    for ii=1:length (unique_IDs)
        ID = unique_IDs(ii);
        tetNo_ID(:,1) = ii;
        indices = tetNo_ID(:,1) == ID;
        unique_tets = unique(tetNo_ID(indices, 2));
        for jj=1: length(unique_tets)
            tet = unique_tets(jj);
            counter=counter + 1;
            indices=find(tetNo_ID(:, 1) == ID & tetNo_ID(:, 2) == tet);
            co_recorded_cell_count = length(indices) - 1;
            output(counter, 1) = ID;
            output(counter, 2) = tet;
            output(counter, 3) = co_recorded_cell_count;
        end
       
    end
    for jj=1: length(tetNo_ID)
        ID = tetNo_ID(jj, 1);
        tet = tetNo_ID(jj, 2);
        tetNo_ID(jj, 3) = output(find(output(:,1) == ID & output(:, 2) == tet) , 3);
    end
    
    corecorded = tetNo_ID(:,3);
%     activeCells = zeros (500,5); 
%     for ii = 1:5
%         corecorded(:,ii) = corecorded(ii)
%          for jj = 1:length(nSpks)
%             if isnan (nSpks (jj,:))
%                 corecorded (jj) = NaN;
%             elseif nSpks(jj,:) > 10
%                 corecorded (jj)= corecorded(jj);
%             else nSpks(jj,:) <= 10;
%                 corecorded (jj)= corecorded(jj) -1;
%             end
%          end
%     end
%      corecorded
end