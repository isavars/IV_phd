function corecorded = corecordedCells(spatData)
%find the number of co-recorded cells per tetrode 
% use tetNo and date to count co-recorded cells per tetrode 
% Also- make a cell # by environment array that shows how many cells were
% corecorded in each environment 
    load ('r889_P18_spatData.mat', 'spatData')
    
    cellID = spatData.cellID; 
%     nSpks = spatData.nSpks; %load spike numbers from file 
%     env = spatData.env; %load environments 
%     nSpks(~(env=='fam' | env=='nov' | env=='sleep')) = NaN;
    
    tetNo_age=zeros(length(cellID), 4); %array to contain cell ID info including co-recorded cells 
    for jj= 1: length(cellID)
        newCellID= string(cellID{jj});
        tetNo_age(jj,1)= double(extractBetween(newCellID,'r','_')); %rat number
        tetNo_age(jj,2)= double(extractBetween(newCellID,'t','c')); %tetrode number
        tetNo_age(jj,3)= double(extractBetween(newCellID,'P','t'));
        if isnan(tetNo_age(jj,3))
            tetNo_age(jj,3)= double(extractBetween(newCellID,'P','_'));%ages
        end
    end
    unique_IDs = unique(tetNo_age(:,1));
    unique_ages = unique(tetNo_age(:,3));
    output = zeros(length(cellID), 4);
    counter=0;
    for ii=1:length (unique_IDs)
        ID = unique_IDs(ii);
        indices = tetNo_age(:,1) == ID;
        unique_ages = unique(tetNo_age(indices, 3));
        for kk=1: length(unique_ages)
            age = unique_ages(kk);
            indices = tetNo_age(:,3) == age;
            unique_tets = unique(tetNo_age(indices, 2));
            for jj=1: length(unique_tets)
                tet = unique_tets(jj);
                counter=counter + 1;
                indices=find(tetNo_age(:, 1) == ID & tetNo_age(:, 2) == tet & tetNo_age(:,3) == age);
                co_recorded_cell_count = length(indices) - 1;
                output(counter, 1) = ID;
                output(counter, 2) = tet;
                output(counter, 3) = age;
                output(counter, 4) = co_recorded_cell_count;
            end
        end
    end
    for jj=1: length(tetNo_age)
        ID = tetNo_age(jj, 1);
        tet = tetNo_age(jj, 2);
        age = tetNo_age(jj, 3);
        tetNo_age(jj, 4) = output(find(output(:,1) == ID & output(:, 2) == tet & output(:, 3) == age), 4);
    end
    
    corecorded = tetNo_age(:,4);
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