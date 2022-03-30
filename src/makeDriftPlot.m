function makeDriftPlot (obj)
    %this makes a raster plot of location of activity (max channels only) on the probe from an nexus object once loaded. 
    % Could adapt to just be from spike times. 
%     TO DO - get the y axis to have channel labels arranged by depth and
%     show blank channels when its not the max channel or get it to display
%     all spikes across all channels 
    dataPath = 'F:\220218i_sleepHP\';

    chanPos              = double(readNPY(fullfile('F:\220218i_sleepHP\','channel_positions.npy'))); 
    chanMapKS            = double(readNPY(fullfile('F:\220218i_sleepHP\','channel_map.npy'))) + 1;
    spkTs                = double(readNPY(fullfile('F:\220218i_sleepHP\','spike_times.npy')));
    clustIDs             = double(readNPY(fullfile('F:\220218i_sleepHP\','spike_clusters.npy'))) + 1;
    clu_info             = tdfread(fullfile('F:\220218i_sleepHP\','cluster_info.tsv'),'tab');

    clu_Depth            = clu_info.depth;%this is the same as the y coords in chaPos but it could be matched onto clu_Ch
    clu_Ch               = clu_info.ch + 1; %then label spkTs with max channels sort by depth and plot all spikes per channel 
    clu_ID                = clu_info.cluster_id + 1;
    clu_Ch               = [clu_ID, clu_Ch, clu_Depth];% make the clust ids match the max channels on the info sheet 
    

%     cluChN = repmat(nan, max(clu_ID),3);
% 
%     % pad clu_Ch with NaN rows for missing clusters 
%     for itcluCh = 1: max(clu_ID)
%         cluChN(itcluCh,1) = itcluCh;
% %          if clu_Ch(cluChN(itcluCh,1),1) == cluChN(itcluCh,1)
% %             cluChN(itcluCh,:) = clu_Ch(cluChN(itcluCh,1),:);
% %          else 
% %          end
%     end
% 
%     for itIDs = 1: length(clu_Ch)
%         if clu_Ch(itIDs,1) == cluChN(itIDs,1)
%             cluChN(itIDs,:) = clu_Ch(itIDs,:);
%         else 
%             cluChN(itIDs+1,:) = clu_Ch(itIDs,:);
%         end
%     end
% 
%     %arrange by depth 
% %     cluChN = sortrows(cluChN,3,'ascend');
% 
% 
%     %make clust_IDs into max channel ids 
%     max_Ch = nan(length(clustIDs),1);
% 
%     for itIDs = 1:length(clustIDs)
%          if cluChN(clustIDs(itIDs),1) == clustIDs(itIDs)
%             max_Ch(itIDs) = cluChN(clustIDs(itIDs),2);
%          else 
%          end
%     end 
%     

% for blank channels in a max channel/good unit only plot add an if
% statement to this bit inside a loop that goes through the channel map top
% to bottom if it has a max channel fill with spike info if not leave
% blank. only issue is spikes with same max channel might be being replaced
% by the second instance 


    spkTs = [];
    chans = [];
    isGood = [obj.cell_ID(:,3), [1:12]'];
    map = [chanMapKS,chanPos];
    map = sortrows(map,3,'ascend');

    %make map index for good clusters 

    [~,map_idx]= ismember (isGood(:,1), map(:,1));
    map_idx = [map_idx,[1:12]'];
    ~ismember (map(:,1)',map(map_idx(:,1)));
    counter =1;
    chCT =1;

    for itMap = map(:,1)'
        if  ismember(itMap,map(map_idx(:,1))) 
            spkTs = [obj.spikeData.spk_Times{1, 1}{map_idx(counter,2), 1}; spkTs]; 
            chans = [repmat(chCT,[length(obj.spikeData.spk_Times{1, 1}{map_idx(counter,2), 1}),1]);chans];%[repmat(map(map_idx(counter,1)),[length(obj.spikeData.spk_Times{1, 1}{map_idx(counter,2), 1}),1]);chans];
            counter = counter +1;
        else
            spkTs = [0; spkTs]; 
            chans = [chCT;chans];%[map(itMap,1);chans];
        end
        chCT = chCT +1;
    end 

%     for itClu = 1:length(obj.spikeData.spk_Times{1, 1})
%         spkTs = [obj.spikeData.spk_Times{1, 1}{itClu, 1}; spkTs]; 
%         chans = [repmat(itClu, [length(obj.spikeData.spk_Times{1, 1}{itClu, 1}),1]);chans];
%     end 
   spkTs = seconds(spkTs);


   s = spikeRasterPlot(spkTs, chans);
end 