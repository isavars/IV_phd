function  getCellProperties_ks()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sampleRate   = 48000;
spikeTimes   = readNPY('F:\r1040_surgery_CA1\spike_times.npy');
spikeTimes   = double(spikeTimes) ./ sampleRate;
% load cluster IDs
clustIDs          = readNPY('F:\r1040_surgery_CA1\spike_clusters.npy') + 1; % 0 based index
 % now we need to remove bad clusters and spike times outside trial
unClustIDs         = unique(clustIDs);
loadFromPhy = 0;
    if loadFromPhy   == 1 
              % phy curated output - only load clusters labelled good
              clu_info       = tdfread('F:\r1040_surgery_CA1\cluster_info.tsv','tab');
              goodLabel    = strcmp(cellstr(clu_info.KSLabel),'good');
              good_clusts    = clu_info.id(goodLabel) + 1;
              clu_Depth      = clu_info.depth(goodLabel);  %
              clu_Ch         = clu_info.ch(goodLabel) + 1; % this is 0 based
              cluLabel       = string(clu_info.group);
              cluLabel       = cluLabel(goodLabel);
    else
            % raw kilosort output - we'll take everything in this case and can
            % filter later if needed
            ks_labels     = tdfread('F:\r1040_surgery_CA1\cluster_KSLabel.tsv','tab'); %
            cluLabel       = string(ks_labels.KSLabel); % this is either 'good' or 'mua'
             %     ind_good       = strcmp(cluLabel,'good');
             good_clusts    = ks_labels.cluster_id + 1; % ID for clusters also 0 based
            % get depth estimate for each cluster based on template ampl
              % distribution across probe
                    templates          = readNPY('F:\r1040_surgery_CA1\templates.npy');
                    Winv               = readNPY('F:\r1040_surgery_CA1\whitening_mat_inv.npy');
                    chanPos            = readNPY('F:\r1040_surgery_CA1\channel_positions.npy');
                    chanMapKS          = double(readNPY('F:\r1040_surgery_CA1\channel_map.npy')) + 1;  %
                    [clu_Depth,clu_Ch] = scanpix.npixUtils.getCluChDepthFromTemplates(templates, Winv, [chanMapKS(:) chanPos(:,2)]);
    end         



% index for 'good' clusters
unGoodClustIDs = unClustIDs(ismember(unClustIDs,good_clusts)); % remove 'mua' or 'noise' clusters from list in case we deal with phy output
% [~,indGood]    = ismember(clustIDs,unGoodClustIDs); % only keep these
% % remove from data
% clustIDs       = clustIDs(~=indGood);
% spikeTimes     = spikeTimes(~=indGood);
% sort clusters, so accumarray output is sorted
[clustIDs, sortInd] = sort(clustIDs);
spikeTimes          = spikeTimes(sortInd);

% reformat into more convenient form
spikeTimesFin  = accumarray(clustIDs,spikeTimes,[max(unGoodClustIDs) 1],@(x) {x}); %max(unGoodClustIDs)

%plots 

%wf plot
% 
[waveforms, channels] = getWaveformsTemp('F:\r1040_surgery_CA1\211122c_surgery.dat',spikeTimesFin(56), 21);
 hold all;

%  figure()
 plot(waveforms{1}(:,:,21)','b'); % for clu56, channel 21 is the max channel! 
%  [waveforms, channels] = getWaveformsTemp('C:\Users\Isabella\Documents\Kilosort_output\Kilosort2\r1040_surgery_CA1\211122c_surgery.dat',spikeTimesFin(56), 16);
%  hold all;
%  figure()
%  plot(waveforms{1}(:,:,29)','b'); % for clu56, channel 21 is the max channel! 


%waveforms = scanpix.npixUtils.extract_waveforms('C:\Users\Isabella\Documents\Kilosort_output\Kilosort2\r1040_surgery_DG\temp_wh.dat',spikeTimesFin{1}, clu_Ch(1))

%AC plot
% %     for it_clu = 1: length(spikeTimesFin)
%              figure()
%             [xc, lags]=spk_crosscorr(spikeTimesFin{38},'AC',0.001,0.1,301,'plot',gca);
% %     end


%     maxRowPerFig = 10;
%     axRowCount = 1;
%     textContent = num2str(clustIDs);
%     for it_clu = 1: length(spikeTimesFin) 
%         if it_clu == 1 || axRowCount > maxRowPerFig
%             axRowCount = 1;
%             hFig = gra_multiplot(maxRowPerFig, 2 );% , 'figborder', [2 1 1 1]
%             axArr = getappdata(hFig, 'axesHandles' ); % makes the axes 
%         end
%         if it_clu == 28
%             figure()
%         else  
%             [xc, lags]=spk_crosscorr(spikeTimesFin{it_clu},'AC',0.002,0.5,300,'plot',axArr(axRowCount,2));
%         end
%         text (axArr(axRowCount,1),-50,23,textContent(it_clu), 'FontSize', 16); 
%         axRowCount = axRowCount + 1;   
% 
%     end


end