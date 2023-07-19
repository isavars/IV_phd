function [cellInfo] = getCellInfo(spatData)

%getCellInfo gets the tetrode number, cell number, age and rat ID from the
%cell ID in spatData being used and produces a matrix with the information
%it also makes a co_recorded_cell_count.

    cellInfo = spatData.cellID;
    dataset = spatData.dataset;
    
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
    
    %this isnt customizable for adults or recordings from a pup on the same 
    %day - needs to inccorporate unique datasets. 

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

end

