function [cellInfo] = getCellInfo()

%getCellInfo gets the tetrode number, cell number, age and rat ID from the
%data file being used and produces a table with the information
%

    cellInfo = loadSpatData('cellID');
    
    tetNo_age=zeros(length(cellInfo), 4); %array to contain cell ID info including co-recorded cells 
    for jj= 1: length(cellInfo)
        newCellInfo= string(cellInfo{jj});
        tetNo_age(jj,1)= double(extractBetween(newCellInfo,'r','_')); %rat number
        tetNo_age(jj,2)= double(extractBetween(newCellInfo,'t','c')); %tetrode number % added extra t for adult trials
        if isnan(tetNo_age(jj,2))
            tetNo_age(jj,2)=  double(extractBetween(newCellInfo,'t t','c'));%for adult trials
        end
%         tetNo_age(jj,3)= double(extractBetween(newCellInfo,'P','t')); % age (40 for adults)40; %
%         if isnan(tetNo_age(jj,3))
%             tetNo_age(jj,3)=  double(extractBetween(newCellInfo,'P','_'));%ages 40; %
%         end
        tetNo_age(jj,3)= double(extractBetween(newCellInfo,'r','_'));% age (this is only a temporary solution that works for this data set FIX IT!)
        if tetNo_age(jj,3) < 900
            tetNo_age(jj,3)= double(extractBetween(newCellInfo,'P','t'));
            if isnan(tetNo_age(jj,3))
                tetNo_age(jj,3)=  double(extractBetween(newCellInfo,'P','_'));%ages 40; %
            end
        else
            tetNo_age(jj,3) = double(extractBetween(newCellInfo,'_','_'));
        end
        tetNo_age(jj,4) = double(extractAfter(newCellInfo,'c')); %co-recorded cells
    end
    
    unique_IDs = unique(tetNo_age(:,1));
    unique_ages = unique(tetNo_age(:,3));
    output = zeros(length(cellInfo), 4);
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
    cellInfo = tetNo_age;

end

