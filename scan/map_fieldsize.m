
function [fs] = map_fieldsize(map, varargin)

% Find field size.
%
%       [fs] = map_fieldsize(map);
%       [fs] = map_fieldsize(map, thr);
%       [fs] = map_fieldsize(map, thr, 'contig');
%
% thr format: 'pc20' = 20% peak rate (default)
%             'hz1'  = 1hz
%             'mean' = Mean rate
%
% 'contig' - field must be a 8-connected contiguous area around peak bin.

% Input Options %
thr_mode = 1;
thr = 0.2;
contig_mode = 0;
if length(varargin)>=1
    thr_str = varargin{1};
    if strcmp(thr_str(1:2), 'pc')
        thr = str2double(thr_str(3:end)) / 100;
        thr_mode = 1;
    elseif strcmp(thr_str(1:2), 'hz')
        thr = str2double(thr_str(3:end));
        thr_mode = 2;
    elseif strcmp(thr_str, 'mean');
        thr_mode = 3;
    else
        error('Threshold format not recognised');
    end
end
if length(varargin)==2
    if strcmp(varargin{2}, 'contig')
        contig_mode = 1;
    end
end

% General map stuff %
if max(max(map))==0
    fs = 0;
    return
end
[a, ind1] = max(map);
[peak_rate, ind2] = max(a);
pkr_r = ind1(ind2);
pkr_c = ind2;
vis_ind = find(map~=-1);
n_bins_vis = length(vis_ind);

% Find threshold values %
if thr_mode==1
    high_rate_ind = find(map>(peak_rate*thr));
elseif thr_mode == 2
    high_rate_ind = find(map>thr);
elseif thr_mode == 3
    m = mean(map(vis_ind));
    high_rate_ind = find(map>m);    
end

% Find Field Sizes %
if contig_mode==0
    
    fs = length(high_rate_ind) / n_bins_vis;
   
elseif contig_mode==1
    
    mapbw = zeros(size(map));
    mapbw(high_rate_ind) = 1;
    mapbw_contig = bwselect(mapbw, pkr_c, pkr_r, 8);
    fs = length(find(mapbw_contig)) / n_bins_vis;
    
end
    