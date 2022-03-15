function [filteredInd,sampleSize,chunkInd] = pos_speedfilter(trialData,varargin)
% Filter data by speed
%
%       [filteredInd,sampleSize,chunkInd] = pos_speedfilter(trialData,varargin)    <with a SCAN structure>
%       [filteredInd,sampleSize,chunkInd] = pos_speedfilter(speedVector,varargin)  <with a vector of speed values. Samp rate assumed to be 50HZ unless specified>
%
% filteredInd is the index into position samples for the speeds that pass.
% sampleSize is the amount of data that passes, in seconds.
% chunkInd labels individual chunks for chunk-filtered data. 
%
% Input options:
%       'min'    - remove speed below a minimum (units are cm/s)
%       'max'    -       "      above a maximum.              
%       'median' - filter speeds to match this median
%       'chunk'  - Only use contiguous chunks of this many seconds

% Parse input %
opt.median = [];
opt.min = 0;
opt.max = 400;
opt.chunk = [];
opt.sample_rate=50;
for ii=1:2:length(varargin)
    opt.(varargin{ii}) = varargin{ii+1};
end
% Allow for just speed vector as main input (assume SR=50, unless specified otherwise) %
if ~isstruct(trialData)
    temp.speed=trialData;
    temp.sample_rate=opt.sample_rate;
    trialData=temp;
end

filteredInd = 1:length(trialData.speed);

%%% Filter by min, max %%%
v = version;
if strcmp(v(1),'8')
    filteredInd = intersect( filteredInd, find(trialData.speed>=opt.min) , 'legacy');  % Between 2010 and 2014, the orientaion of the output
    filteredInd = intersect( filteredInd, find(trialData.speed<=opt.max) , 'legacy');  % of intersect has changed. Don't have time to fix properly now.
else
    filteredInd = intersect( filteredInd, find(trialData.speed>=opt.min) );  
    filteredInd = intersect( filteredInd, find(trialData.speed<=opt.max) );  
end

%%% Filter to match a median speed %%%
if ~isempty(opt.median)
    speeds = trialData.speed(filteredInd);
    % This code is adapted from medianAdjust (AJ) %
    [orderedSpeeds,orderedSpeedIndex] = sort(speeds);
    % Find the speed and pos sample of the speed which is closest to the grand median for this rat
    absDiff = abs(speeds-opt.median);
    trialGrandMedianIndex = find(absDiff==min(absDiff),1);
    % Find index into orderedSpeedIndex of the grand median speed (in trial)
    orderedTrialGrandMedianIndex = find(orderedSpeedIndex==trialGrandMedianIndex,1);
    % Find which chunk of data above and below trial grand median speed is the smaller
    % and use that to select trial grand median +/- the smaller chunk
    belowTrialGrandMedianSize = orderedTrialGrandMedianIndex-1;
    aboveTrialGrandMedianSize = length(speeds)-orderedTrialGrandMedianIndex;
    % Find allowed data - min and max speed so that distribution between these matches median %
    allowedDataSize = min([belowTrialGrandMedianSize,aboveTrialGrandMedianSize]);
    minAllowedSpeed = speeds(orderedSpeedIndex(orderedTrialGrandMedianIndex-allowedDataSize));
    maxAllowedSpeed = speeds(orderedSpeedIndex(orderedTrialGrandMedianIndex+allowedDataSize));
    % Find index of speeds within range. Note index is into original, unfiltered speed data %
    medianFilterInd = find( trialData.speed>=minAllowedSpeed & trialData.speed<=maxAllowedSpeed );
    % Combine with previous filters %
    filteredInd=intersect(filteredInd,medianFilterInd);
end

%%% Take chunks longer than X sec threshold %%%
% This is copied from posChunks (AJ) %
if ~isempty(opt.chunk)
    contiguousChunks = diff(filteredInd);
    endsOfContiguousChunks = [find(contiguousChunks>1), length(filteredInd)];
    lengthOfContiguousChunks = [endsOfContiguousChunks(1), diff(endsOfContiguousChunks)]; %CB potentially miss first continguous chunk?
    allowedContiguousChunks = find(lengthOfContiguousChunks>=(opt.chunk*trialData.sample_rate));
    endsOfAllowedContiguousChunks = endsOfContiguousChunks(allowedContiguousChunks);
    lengthOfAllowedContiguousChunks = lengthOfContiguousChunks(allowedContiguousChunks);
    chunkSamples = zeros(max(lengthOfContiguousChunks),length(lengthOfAllowedContiguousChunks)); % The indices of each allowed chunk will be a column.
    chunkInd = zeros(max(lengthOfContiguousChunks),length(lengthOfAllowedContiguousChunks));
    for kk = 1:length(lengthOfAllowedContiguousChunks)
        temp = filteredInd(endsOfAllowedContiguousChunks(kk)-lengthOfAllowedContiguousChunks(kk)+1):filteredInd(endsOfAllowedContiguousChunks(kk));
        chunkSamples(1:length(temp),kk) = temp';
        chunkInd(1:length(temp),kk) = kk;
    end
    filteredInd = chunkSamples(chunkSamples~=0); % Convert to linear index
    chunkInd = chunkInd(chunkSamples~=0);
else
    chunkInd = ones(size(filteredInd));
end

sampleSize = length(filteredInd)/trialData.sample_rate;





