%% histInfo - labels each tetrode* with the most likely layer it was in based on histology
% reads histology table info and cellID info from spatData and outputs a label: GCL, HL, CA3, GCL HL, 
% HL CA3 or CA3 HL and a label for probe types for each cell based on it's tetrode to be used by any function. 
% *labelling per tetrode for now as part of the SfN data push 

function [cell_layer]= histInfo()

    load ('histology.mat', 'histology')

    %get tet data from cell info
    cellInfo = getCellInfo();
    rat = [1099;222;1241;1242];%cellInfo(:,1);%
    tet = [2;3;4;7]; %cellInfo(:,2);%

    %loop through animal ids to see if they have histology and add the
    %histologylabel per tetrode to each cell (lengthwise on spat data)
    %column one of cell layer contains the hitology label, column 2
    %contains the cetainty raiting and column 3 contains the probe type. 
    
    cell_layer = cell(length(rat),3);

    for it_ID = 1:length(rat)
        if ismember(rat(it_ID), histology.rat_ID)
            tetrode_number = tet(it_ID);
            rat_row = find(histology.rat_ID == rat(it_ID));
            if tetrode_number == 1 
                cell_layer{it_ID,1} =  histology.channel_group_1(rat_row);
                if contains(cell_layer{itID,1}," ") == 1
                    cell_layer{it_ID,2} = 2; %certainty score of 2 if there are two possible layers
                else
                    cell_layer{it_ID,2} = 1;  %certainty score of 1 if there is one possible layer
                end
            elseif tetrode_number == 2
                cell_layer{it_ID,1} = histology.channel_group_1(rat_row);
                if contains(cell_layer{it_ID,1}," ") == 1
                    cell_layer{it_ID,2} =  2;
                else
                    cell_layer{it_ID,2} = 1;
                end
            elseif tetrode_number == 3
                cell_layer{it_ID,1} = histology.channel_group_2(rat_row);
                if contains(cell_layer{it_ID,1}," ") == 1
                    cell_layer{it_ID,2} =  2;
                else
                    cell_layer{it_ID,2} = 1;
                end
            elseif tetrode_number == 4 
                cell_layer{it_ID,1} = histology.channel_group_2(rat_row);
                if contains(cell_layer{it_ID,1}," ") == 1
                    cell_layer{it_ID,2} =  2;
                else
                    cell_layer{it_ID,2} = 1;
                end
            elseif tetrode_number == 5
                cell_layer{it_ID,1} = histology.channel_group_3(rat_row);
                if contains(cell_layer{it_ID,1}," ") == 1
                    cell_layer{it_ID,2} =  2;
                else
                    cell_layer{it_ID,2} = 1;
                end
            elseif tetrode_number == 6 
                cell_layer{it_ID,1} = histology.channel_group_3(rat_row);
                if contains(cell_layer{it_ID,1}," ") == 1
                    cell_layer{it_ID,2} =  2;
                else
                    cell_layer{it_ID,2} = 1;
                end
            elseif tetrode_number == 7
                cell_layer{it_ID,1} = histology.channel_group_4(rat_row);
                if contains(cell_layer{it_ID,1}," ") == 1
                    cell_layer{it_ID,2} =  2;
                else
                    cell_layer{it_ID,2} = 1;
                end
            elseif tetrode_number == 8
                cell_layer{it_ID,1} = histology.channel_group_4(rat_row);
                if contains(cell_layer{it_ID,1}," ") == 1
                    cell_layer{it_ID,2} =  2;
                else
                    cell_layer{it_ID,2} = 1;
                end
            end
        else
            cell_layer{it_ID,1} = nan;
            

        end
        cell_layer{it_ID,3} = histology.probe_type(rat_row);

    end


    
end 