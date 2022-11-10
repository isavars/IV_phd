function [layer] = makeDS2(eegs,mapping)
%I want this to: 
    % 1. load eegs from selected trials 
    % 2. threshold for DS2 spike and select strongest ampliudes - combine this with actual
    % depth from probe measurements 
    % 3. make a figure showing the DS2 spike organized by depth on the
    % probe. - just for reference - can plot depth vs amplitude of spike to
    % see if it goes from positive to negative - this could determine the
    % placement. 
    % 4. determine if pointing up or down (pre or post granule cell layer)-
    % bear in mind signal recorded on silicone probes will be inverted?
    % check this. 
    % 5. output a label per channel (or per tetrode -SfN) that indicates 
    % likelihood of being in the GCL or HL - at or below DS2 reversal 
%     notes from Senzai and Buzsaki 2017 - Dentate spike 2 (DS2) was first 
    % detected as the period when the difference of the band pass filtered 
    % (2-50Hz) LFP recorded from a hilus location and molecular layer site 
    % exceeded 1.14mV. If the mean LFP value at the molecular layer during the
    % DS2 was lower than the baseline value (−36 ms to −16 ms) by > 0.19mV, 
    % the DS2 passed the criteria and included for later analysis. Peak time of 
    % DS2 was defined as the time when the dentate hilus wide-band LFP showed 
    % maximum value.

    eeg_output = [];
    for it_eegs = 1:32
        if it_eegs == 1
            eeg_struct  = load_eeg(strcat('220218i_sleepHP_LP.eeg'));
        else
            eeg_struct = load_eeg(strcat('220218i_sleepHP_LP.eeg',num2str(it_eegs)));
        end
        eeg_output = [eeg_struct.eeg;eeg_output];
    end
    
    eeg_output = transpose(reshape(eeg_output,[300000,32])); %find out what that number is so it's not hc
    
    %convert to mV 
    
    voltages = eeg_output.*(1.5 ./ 128);
    voltages = voltages.';
    
    %make give each channel a depth label - using depth map tables 

    load ('mapOpps.mat', 'singleMapOpps', 'multiMapOpps')
    
        %find histology labels this makes a list of labels of row length
        %getSpatData so per cell. Might get the label like this in the future
        %but simplfying for now 
    
    %     cell_layer= histInfo();
    %     probe_type = cell_layer{:,3};
    % 
    %     %chose what map to use per cell based on probe type 
    %     for itCell = 1:length(probe_type)
    %         if strcmp(probe_type (itCell), "1x32")
    %             % Re-ordering based on CM32 1x32 plugged in with NN and Omnetics on opposite sides
    %             voltages = [voltages(8,:);voltages(9,:);voltages(16,:);voltages(1,:);voltages(7,:);voltages(19,:);voltages(30,:);voltages(10,:);voltages(15,:);voltages(2,:);voltages(25,:);voltages(20,:);voltages(29,:);voltages(24,:);voltages(14,:);voltages(3,:);voltages(26,:);voltages(21,:);voltages(28,:);voltages(23,:);voltages(13,:);voltages(4,:);voltages(32,:);voltages(22,:);voltages(27,:);voltages(17,:);voltages(12,:);voltages(5,:);voltages(31,:);voltages(11,:);voltages(6,:);voltages(18,:)];
    %         elseif strcmp(probe_type (itCell), "4x8")
    %              % Re-ordering based on CM32 Buzsaki 4x8 plugged in with NN and
    %              % Omnetics on opposite sides - figure out which one is right!!
    %             voltages = [voltages(6,:);voltages(2,:);voltages(5,:);voltages(29,:);voltages(27,:);voltages(3,:);voltages(4,:);voltages(28,:);voltages(30,:);voltages(8,:);voltages(1,:);voltages(7,:);voltages(31,:);voltages(25,:);voltages(32,:);voltages(26,:);voltages(9,:);voltages(19,:);voltages(10,:);voltages(16,:);voltages(24,:);voltages(18,:);voltages(23,:);voltages(17,:);voltages(15,:);voltages(11,:);voltages(20,:);voltages(12,:);voltages(14,:);voltages(22,:);voltages(21,:);voltages(13,:)];
    %         end
    %     end 
    
    
            if strcmp(probe_used, "single")
                % Re-ordering based on CM32 1x32 plugged in with NN and Omnetics on opposite sides
                disp('single')
    %             voltages = [voltages(8,:);voltages(9,:);voltages(16,:);voltages(1,:);voltages(7,:);voltages(19,:);voltages(30,:);voltages(10,:);voltages(15,:);voltages(2,:);voltages(25,:);voltages(20,:);voltages(29,:);voltages(24,:);voltages(14,:);voltages(3,:);voltages(26,:);voltages(21,:);voltages(28,:);voltages(23,:);voltages(13,:);voltages(4,:);voltages(32,:);voltages(22,:);voltages(27,:);voltages(17,:);voltages(12,:);voltages(5,:);voltages(31,:);voltages(11,:);voltages(6,:);voltages(18,:)];
                %assign channels to ycoords for single using singleMapOpps
            elseif strcmp(probe_type, "multi")
                 % Re-ordering based on CM32 Buzsaki 4x8 plugged in with NN and
                 % Omnetics on opposite sides - figure out which one is right!!
                disp('multi')
    %              voltages = [voltages(6,:);voltages(2,:);voltages(5,:);voltages(29,:);voltages(27,:);voltages(3,:);voltages(4,:);voltages(28,:);voltages(30,:);voltages(8,:);voltages(1,:);voltages(7,:);voltages(31,:);voltages(25,:);voltages(32,:);voltages(26,:);voltages(9,:);voltages(19,:);voltages(10,:);voltages(16,:);voltages(24,:);voltages(18,:);voltages(23,:);voltages(17,:);voltages(15,:);voltages(11,:);voltages(20,:);voltages(12,:);voltages(14,:);voltages(22,:);voltages(21,:);voltages(13,:)];
                %assign channels to ycoords for multiMapOpps 
                multiMapOpps = sortrows(multiMapOpps,4);%check what colun has the correct orentation
                ycoords = multiMappOpps.ycoords;
            end
    
        %give each channel a depth label using coordinates in mapOpps table 

    
        
        %absolute_depth =  

        %basically run something like extract spikes for all 32 channels
        %using the threshold value 1.14 mV - when there's a crossing fill
        %in a matrix with 50 ms of data peak should be in the middle - exclude 
        % anything where the points at either side dont decrease back to
        % baseline 
        % add 1-32 voltages at the same timestamp - this is something that 
        % you can plot as a sense check by depth using the labels. - take each
        %of these voltages (the peak and the surrounding points across the
        %32 channels and filter for an inversion in amplitude - this should
        %confirm if any of the crossings are a true DS2 - take the peak DS2
        %value from the whole recording and use this to find the channels
        %with the inverison point - label channels at inversion point GCL
        %and those below HL (maybe grade them with a distance from inverison metric) 


    
    
    %make big plot with all the eegs arranged by depth for single and multi
    %shank probes 
    
    %single shank plot 
    
        hFig = gra_multiplot(32,1);
        axArr = getappdata(hFig, 'axesHandles');
    
    
    
    
    % plot(eeg_output(1,[500:1000]))
    
%         all_thresh_cross = [];
%         thresh = 1.2;
%     
%         absolute_eeg_output = abs(eeg_output);
%         thresh_cross = absolute_eeg_output > thresh;
%         
%     
%     for it_eeg = 1:300000
%             if (eeg_output(1,it_eeg) > thresh && eeg_output(1,it_eeg + 5) > thresh && eeg_output(1,it_eeg + 10) > thresh) || (eeg_output(1,it_eeg) < -thresh && eeg_output(1,it_eeg + 5) < -thresh && eeg_output(1,it_eeg + 10) < -thresh)
%             thresh_cross = eeg_output(1,[it_eeg-14:it_eeg+17]);
%             all_thresh_cross = [all_thresh_cross;thresh_cross];
%         end
%         
%     end 
    
    
    
    % 
    % all_thresh_cross = all_thresh_cross(1:20:end,:);
    % whos
    % 
    % for it_th = 1:length(all_thresh_cross)
    %     figure()
    %     plot (all_thresh_cross(it_th,:))
    % 
    % end