function [cellStruct]=spk_filter(cellStruct,filtInd,sampleRate)
% Filter spikes, using index obtained from position or EEG data.
%
%       [cellStructFilt]=spk_filter(cellStruct,filtInd,sampleRate)
%
% filtInd should be numerical. sampleRate refers to the data on which the
% filtInd is based (i.e. pos or EEG)

% disp('(24/09/09) Warning!: I think there is a bug in SPK_FILTER, line 13 - using INTERSECT means only first spike in pos bin is used. Use ISMEMBER instead?');
for ii=1:length(cellStruct)
    % Get filter index %
    spkSampInd=ceil(cellStruct(ii).st .* sampleRate);
    cellFiltInd=find(ismember(spkSampInd,filtInd));
    % Filter spike times. Need to adjust times to compensate for (pos or EEG) samples filtered out %
    spkTemp=[];
    for kk=1:length(cellFiltInd)
        oldTime=cellStruct(ii).st( cellFiltInd(kk) );
        oldSamp=spkSampInd( cellFiltInd(kk) );              % The sample that the spike fired in originally.
        lostSamps = oldSamp - find( oldSamp==filtInd, 1 );  % Get the difference between above, and the difference between its original position, and its new position in the filtered data.
        lostTime = lostSamps/sampleRate;                    % Convert the difference between samples to a time difference.
        spkTemp = [spkTemp oldTime-lostTime];
    end
    cellStruct(ii).st = spkTemp;
    % Filter waveforms %
    if ~isempty(cellStruct(ii).wf_amps)
        cellStruct(ii).wf_amps = cellStruct(ii).wf_amps(cellFiltInd,:);
    end
    if ~isempty(cellStruct(ii).wf_all)
        cellStruct(ii).wf_all = cellStruct(ii).wf_all(cellFiltInd,:,:);
    end
end
    
    


