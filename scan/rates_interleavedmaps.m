function [rtnOdd rtnEven]=rates_interleavedmaps(data,win,varargin)
% Interleaved (odd-even) time window maps.
%
%   [mapsOdd mapsEven] = rates_interleavedmaps(data,windowDuration);
%   [mapsOdd mapsEven] = rates_interleavedmaps(data,windowDuration,'paramsProp','paramsValue', .. );

for ii=1:length(data.trials)
    sr = data.trials(ii).sample_rate;
    nWin=( floor( data.trials(ii).dur/(win*2) ) )*2; % So both IL maps have the same number of window chunks.
    winEnd = win:win:win*nWin;
    winStart = [1/sr winEnd(1:end-1)];
    % Generate direct index filters for each map %
    mapOddInd=[];
    for jj=1:2:nWin
        mapOddInd = [mapOddInd ceil(sr*winStart(jj)):floor(sr*winEnd(jj))];
    end
    mapEvenInd=[];
    for jj=2:2:nWin
        mapEvenInd = [mapEvenInd ceil(sr*winStart(jj)):floor(sr*winEnd(jj))];
    end
    % Make rate maps %
    prms = rates_params('filt_index',mapOddInd,varargin{1:end});
    rtnOdd = rates_main(data,prms);
    prms = rates_params('filt_index',mapEvenInd,varargin{1:end});
    rtnEven = rates_main(data,prms);    
end  
        