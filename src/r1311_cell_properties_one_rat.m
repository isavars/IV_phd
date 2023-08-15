%activity tracked across days fro same rat (r1311)
load ("r1311_all_shuff90_spatData.mat", 'spatData')
load('electrode_positions_and_ds2_for_all_rats_5.mat', 'elePos')

cellInfo = getCellInfo(spatData);

%chose trial with most spikes to use wfs, spike times, mean firing rate, max wf channel from it 
for itSp = 1: length (spatData.nSpks) 
    [~, maxSpksPos] = max(spatData.nSpks(itSp,1:5)); %max for all wake trials
    max_wf_chan(itSp,:) = spatData.max_wf_channel(itSp,maxSpksPos); 
end

%I want to lookk at changes in different properties over age - proportions
%cant be looked at because these could be changing due to cell loss on the
%probe - only indivicual cell properties - mean rate, burst index, wf width
%etc and if a cell remaps/is silent etc 
%what im trying to do here is see if theres a difference between shank 1
%which is likely only granule cells and the other shanks (shank 3 is the
%most likely to be hilus 4 could have CA3 and 2 is GCL hilus boarder 

%filter data by shank - make tetshankchan 

[tetShankChan, shank_channels] = makeTetShankChan(spatData,cellInfo, max_wf_chan,elePos);

%i want a plot of mean rate per shank over age




function [tetShankChan, shank_channels] = makeTetShankChan(spatData,cellInfo, max_wf_chan,elePos) 
    %make tetrode labels from CellIDs in spatData - needs to be adaptable
    %for single shank probe data also!!   

        tet = cellInfo(:,2);
        tetShankChan = zeros(length(tet),3);
        shank_channels = zeros(length(tet),8); %this is going to contain the channels that are on each shank to be used by electrode position measure
        % this is kept the same per probe type - in the single shank probes
        % cells recorded on the "same shank" will be from the 8 closest 
        % contacts to the 'tetrode' the channel is on 
        for it_tet = 1: length(tet)            
            % create channel index per tetrode 'tet_index' which is a tetrode by
            % channel identity array for each type of probe 
            if elePos.probe_type(spatData.animal(it_tet) == elePos.rat_ID) == 4 %this creates an index for the row in elepos containing the rat that the cell came from in tet. 
                tet_index = zeros(8, 4); % create an 8x4 matrix of zeros
                tet_index(1:2:end, :) = repmat(1:8:25, 4, 1).' + repmat(0:2:6, 4, 1); % fill odd-numbered rows with odd numbers
                tet_index(2:2:end, :) = repmat(2:8:26, 4, 1).' + repmat(0:2:6, 4, 1); % fill even-numbered rows with even numbers
                overlap = repmat (0:8:24,4,1).' + repmat (5:8,4,1);%overlap tetrodes are the bottom 4 contacts which are closer to each other 
                tet_index = [tet_index; overlap]; % + 1 overlap per octrode    
                %make tetShankChan for multishank probes 
                tetShankChan(it_tet,1) = tet(it_tet);
                tetShankChan(it_tet,3) = tet_index(tet(it_tet),max_wf_chan(it_tet));
                if tet(it_tet) == 1 || tet(it_tet) == 2 || tet(it_tet) == 9
                    tetShankChan(it_tet,2) = 1; %adds shank label 
                    shank_channels(it_tet,:) = reshape(tet_index(1:2, :), [],1)'; %says all the channels on that shank
                elseif tet(it_tet) == 3 || tet(it_tet) == 4 || tet(it_tet) == 10
                    tetShankChan(it_tet,2) = 2;
                    shank_channels(it_tet,:) = reshape(tet_index(3:4, :), [],1)';
                elseif tet(it_tet) == 5 || tet(it_tet) == 6 || tet(it_tet) == 11
                    tetShankChan(it_tet,2) = 3;
                    shank_channels(it_tet,:) = reshape(tet_index(5:6, :), [],1)';
                elseif tet(it_tet) == 7 || tet(it_tet) == 8 || tet(it_tet) == 12
                    tetShankChan(it_tet,2) = 4;
                    shank_channels(it_tet,:) = reshape(tet_index(7:8, :), [],1)';
                end
            elseif elePos.probe_type(spatData.animal(it_tet) == elePos.rat_ID) == 1
                tet_index = 1:32;
                tet_index = reshape(tet_index,4,8).';
                overlap_rows = repmat (0:8:24,4,1).' + repmat (3:6,4,1);
                % Add the overlapping rows to the matrix
                tet_index = [tet_index; overlap_rows];
                % make shank channel index for making shank channels that
                % represent 8 closest contacts to recorded tetrodes
    
                if tet(it_tet) == 1 || tet(it_tet) == 9
                    shank_channels(it_tet,:) = 1:8;
                elseif tet(it_tet) == 8
                    shank_channels(it_tet,:) = 25:32;
                else
                    row = tet_index(tet(it_tet),:);
                    shank_channels(it_tet,:) = [row(1)-2, row(1)-1, row, row(end)+1, row(end)+2];
                end

                %make tetShankChan for single shank - in this case shank
                %represents the 4 closest contacts to the 'tetrode' 
                tetShankChan(it_tet,1) = tet(it_tet);
                tetShankChan(it_tet,3) = tet_index(tet(it_tet),max_wf_chan(it_tet));
                if tet(it_tet) == 1 || tet(it_tet) == 2 || tet(it_tet) == 9
                    tetShankChan(it_tet,2) = 1;
                    %shank_channels(it_tet,:) = reshape(tet_index(1:2, :), [],1)';
                elseif tet(it_tet) == 3 || tet(it_tet) == 4 || tet(it_tet) == 10
                    tetShankChan(it_tet,2) = 2;
                    %shank_channels(it_tet,:) = reshape(tet_index(3:4, :), [],1)';
                elseif tet(it_tet) == 5 || tet(it_tet) == 6 || tet(it_tet) == 11
                    tetShankChan(it_tet,2) = 3;
                    %shank_channels(it_tet,:) = reshape(tet_index(5:6, :), [],1)';
                elseif tet(it_tet) == 7 || tet(it_tet) == 8 || tet(it_tet) == 12
                    tetShankChan(it_tet,2) = 4;
                    %shank_channels(it_tet,:) = reshape(tet_index(7:8, :), [],1)';
                end
            end 

        end
end 