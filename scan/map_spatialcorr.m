
function [r] = map_spatialcorr(rates1, rates2, varargin)

% Correlate two place maps
%
%   r = map_spatialcorr(map1, map2);
%   r = map_spatialcorr(map1, map2, 'ignoreMutualZero');
%
% Use only mutually visited bins. Background must be set to NaN.
% If called with 'ignoreMutualZero' switch, correlation will exclude postion bins
% with a rate of exactly 0Hz in both trials.

if isempty(rates1) || isempty(rates2)
    r = NaN;  return
end
if max(rates1(:))==0 || max(rates2(:))==0
    r = NaN;  return
end

rates1 = double(rates1);   
rates2 = double(rates2);

visInd = ~isnan(rates1) & ~isnan(rates2);   % Find unvisited bins in either trial.

if ~isempty(varargin) 
    if strcmp(varargin{1},'ignoreMutualZero')
        zeroInd = rates1==0 & rates2==0;
        visInd = visInd & ~zeroInd;
    else
        error('Input argument not recognized');
    end
end

r=corr2(rates1(visInd),rates2(visInd));