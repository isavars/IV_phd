function [cellInfo] = getCellInfo(spatData)

%getCellInfo gets the tetrode number, cell number, age and rat ID from the
%cell ID in spatData being used and produces a matrix with the information
%it also makes a co_recorded_cell_count.

    cellInfo = spatData.cellID;
    dataset = spatData.dataset;

    %get values for sleep spikes threshold filtering 
    nSpks = spatData.nSpks;
    %make sleep trial index 
    sleep_idx = zeros(size(spatData,1),1);
    for itCl = 1: height(spatData)
        sleep_trials = strcmp(string(spatData.env(itCl,:)),'sleep');
        sleep_idx_temp = find(sleep_trials);
        if size(sleep_idx_temp,2) > 1 
            sleep_idx(itCl) = sleep_idx_temp(2); %dealing with trial with more than one sleep
        else 
            sleep_idx(itCl) = sleep_idx_temp;
        end
    end

    %make matrix from info extraced from cell ID 
    tetNo_age=zeros(length(cellInfo), 5); %array to contain cell ID info including co-recorded cells 
    for jj= 1: length(cellInfo)
        newCellInfo= string(cellInfo{jj});
        tetNo_age(jj,1)= double(extractBetween(newCellInfo,'r','_')); %rat number
        tetNo_age(jj,2)= double(extractBetween(newCellInfo,'t','c')); %tetrode number % added extra t for adult trials
%         if isnan(tetNo_age(jj,2))
%             tetNo_age(jj,2)=  double(extractBetween(newCellInfo,'t t','c'));%for adult trials
%         end
        tetNo_age(jj,3)= double(extractBetween(newCellInfo,'P','t')); % age (40 for adults)
%         if isnan(tetNo_age(jj,3))
%             tetNo_age(jj,3)=  double(extractBetween(newCellInfo,'P','_'));%ages 40; % is this for age 40 or for trials that are more than one on the same day? 
%         end
        tetNo_age(jj,4) = double(extractAfter(newCellInfo,'c')); %co-recorded cells 
        tetNo_age(jj,5) = double(extractBetween(newCellInfo,'_','_'));

    end
     

    unique_IDs = unique(tetNo_age(:,1));
    output = zeros(length(cellInfo), 4);
    counter=0;
    for ii=1:length (unique_IDs)
        ID = unique_IDs(ii);
        indices = tetNo_age(:,1) == ID;
        unique_ages = unique(tetNo_age(indices, 3));
        for kk=1: length(unique_ages)
            age = unique_ages(kk);
            indices = tetNo_age(:,3) == age;
            unique_dates = unique(tetNo_age(indices,5));
            for ll = 1:length(unique_dates)
                dates = unique_dates(ll); 
                indices = tetNo_age(:,5) == dates;
                unique_tets = unique(tetNo_age(indices, 2));
                for jj=1: length(unique_tets)
                    tet = unique_tets(jj);
                    counter=counter + 1;
                    indices=find(tetNo_age(:, 1) == ID & tetNo_age(:, 2) == tet & tetNo_age(:,3) == age & tetNo_age(:,5) == dates);
                    co_recorded_cell_count = length(indices) - 1; 
                    output(counter, 1) = ID;
                    output(counter, 2) = tet;
                    output(counter, 3) = age;
                    output(counter, 4) = co_recorded_cell_count; 
                    output(counter, 5) = dates;
                end
            end
        end
    end
    
    for jj=1: length(tetNo_age)
        ID = tetNo_age(jj, 1);
        tet = tetNo_age(jj, 2);
        age = tetNo_age(jj, 3);
        dates = tetNo_age(jj, 5);
        tetNo_age(jj, 4) = output(find(output(:,1) == ID & output(:, 2) == tet & output(:, 3) == age & output(:, 5) == dates), 4); %co-recorded-cells 
    end

    cellInfo = tetNo_age;

    % New loop to update values based on nSpks condition
    for jj=1:length(cellInfo)
        ID = cellInfo(jj, 1);
        tet = cellInfo(jj, 2);
        age = cellInfo(jj, 3);
        dates = cellInfo(jj, 5);      
        % Find indices of rows with the same ID, tet, age, and date
        matchingIndices = find(cellInfo(:,1) == ID & cellInfo(:,2) == tet & cellInfo(:,3) == age & cellInfo(:,5) == dates);
        % Check your condition on nSpks. As an example, let's assume you're checking if nSpks(jj) > threshold
        threshold = 50; % This is just an example threshold. You can replace this with your actual condition.
        if nSpks(jj, sleep_idx(jj)) < threshold && cellInfo(jj, 4) ~=0
            % Update all rows (column 4) with matching indices
            cellInfo(matchingIndices, 4) = cellInfo(matchingIndices, 4) - 1; % Replace 'someValue' with the value you want to set.
        end
    end

end

