function [bits_per_spike,bits_per_sec] = map_skaggsinfo(rates, times)

% Find spatial information (skaggs info) of place field
%
%   [bits_per_spike bits_per_second] = map_skaggsinfo(rate_map, pos_map);
%
% Returns Skaggs et al's estimate of spatial information in bits per second:
%
%   I = sum_x p(x) r(x) log(r(x)/r)  
%
% then divide by mean rate over bins to get bits per spike.
%
% Binning could be over any single spatial variable (e.g. location, direction, speed).

if max(max(rates)) <= 0
    bits_per_spike = NaN;
    bits_per_sec = NaN;
    return
end

if (size( rates ) ~= size(times))
    error('Size mismatch between spikes and position inputs to skaggs_info');
end

rates(isnan(rates)) = 0;   % remove background
times(isnan(times)) = 0;

if(size( rates ) > 1)
    % turn arrays into column vectors
    rates = reshape( rates, prod(size(rates)), 1);
    times = reshape( times, prod(size(times)), 1);
end

duration = sum(times);
mean_rate = sum(rates.*times)./duration;

p_x = times./duration;
p_r = rates./mean_rate;                   
dum = p_x.*rates;        
ind = find( dum > 0 );   
bits_per_sec = sum(dum(ind).*log2(p_r(ind)));   % sum( p_pos .* rates .* log2(p_rates) )
bits_per_spike = bits_per_sec/mean_rate;

