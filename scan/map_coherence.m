
function [rtn] = map_coherence(map, varargin)

% Spatial coherence of a rate map.
%
%       [rtn] = map_coherence(map);
%       [rtn] = map_coherence(map, N);
%
% Correlates rate of bin, and mean rate of the surrounding
% bins, over all visited bins. Maps should be non-smoothed.
%
% N is the number of surrounding bins to consider, expressed
% the length of the side of the square of the bin + surrounding.
% Default = 3 (i.e. 8 surrounding bins).
%
% e.g. Kubie, Muller & Bostock (1990).

if max(max(map))==0
    rtn = NaN;
    return
end
if length(varargin)==1
    boxcar = varargin{1};
else
    boxcar = 3;
end

if any(size(map)==1)
    %% Vector maps, for heading correlates %%
    dim=1;
    % Get important indices from neighbourhood size %
    bin_ind = ceil((boxcar^dim) / 2);
    nhood_ind = setdiff(1:boxcar^dim, bin_ind);
    n_pad = floor(boxcar/2);
    % Columnise map - pad first, so that all bins are included in colmap, and
    % therefore indices match up %
    pad = [map(end-n_pad+1:end); map; map(1:n_pad)];
    colmap = im2col(pad, [boxcar 1], 'sliding');
    % Find rates of visited bins, and mean of 8 neighbourhood bins (
    bins = colmap(bin_ind,:);
    nhood = colmap(nhood_ind,:);
    nhood = nanmean(nhood);
    % Correlate %
    rtn = corr2(bins, nhood);
else
    %% 2-D Place Maps %%
    % Get important indices from neighbourhood size %
    bin_ind = ceil((boxcar^2) / 2);
    nhood_ind = setdiff(1:boxcar^2, bin_ind);
    n_pad = floor(boxcar/2);
    % Columnise map - pad first, so that all bins are included in colmap, and
    % therefore indices match up %
    s = size(map);
    pad = repmat(NaN, s+(n_pad*2));
    pad((1+n_pad):(s(1)+n_pad), (1+n_pad):(s(2)+n_pad)) = map;
    colmap = im2col(pad, [boxcar boxcar], 'sliding');
    % Find rates of visited bins, and mean of 8 neighbourhood bins (if visited)
    colmap = colmap(:,~isnan(map));
    bins = colmap(bin_ind,:);
    nhood = colmap(nhood_ind,:);
    nhood(nhood==-1) = NaN;
    nhood = nanmean(nhood);
    % Correlate %
    rtn = corr2(bins, nhood);
end



