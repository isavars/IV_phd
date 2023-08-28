function []=spk_psth(spikeTrain, stimTimes, window, binSize, varargin)
% Peri-stimulus time histogram of spikes.
%
%       spk_psth(spikeTrain, stimTimes, window, binSize)
%       spk_psth(spikeTrain, stimTimes, window, binSize, axisHandle)
%
% spikeTrain and stimTimes should be specified in seconds (from trial initiation).
%
% Specify window (i.e. start and end for histogram, relative to stimulus onset)
% and binSize in ms.
%
% if axesHandle argumnet is given, will plot into specified axis.


window=window./1000;    binSize=binSize/1000;

binEdgeTimes = window(1) : binSize : window(2);

% Construct index for filtering particular stimulation episodes %
filtWindow=[window(1) 0; 2 window(2)];
histFiltIndex=zeros(size(binEdgeTimes));
for ii=1:size(filtWindow,1)
    histFiltIndex( binEdgeTimes>=filtWindow(ii,1) & binEdgeTimes<filtWindow(ii,2) ) = 1;
end
histFiltIndex=logical(histFiltIndex);

% Get histogram for each episode %
histByStim = nan(length(stimTimes),length(binEdgeTimes));
episodesUsed=0;
for ii=1:length(stimTimes)
    binTemp = binEdgeTimes + stimTimes(ii);
    temp = histc(spikeTrain,binTemp);
    if sum( temp(histFiltIndex') ) >= 5
        histByStim(ii,:) = temp;
        episodesUsed=episodesUsed+1;
    end    
end

% Get mean histogram %
meanHist=nansum(histByStim)./episodesUsed;

if isempty(varargin)
    figure;
    bar(binEdgeTimes,meanHist,'histc');
    shading(gca,'flat');
else
    hAx=varargin{1};
    hFig=get(hAx,'parent');
    set(hFig,'renderer','painters');
    bar(hAx,binEdgeTimes,meanHist,'histc');
    shading(hAx,'flat');
    hold(hAx,'on');
    plot(hAx, window, [0 0], 'k-');
end