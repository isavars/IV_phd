function [rtn]=map_fft2(map)
% This function does not work - I got stuck thinking about how to pass the frequency output.
% Maybe not worth the bother?

% 2-dimensional Fourier Transform power spectrum of a rate map.
%
%   PS = map_fft2(map);
%   [PS spatial_frequencies] = map_fft2(map,bins_per_metre,crop_factor);
%
% Specifying a 'crop_factor' argument will crop the power spectrum to
% the dimension map_size*crop_factor.


map=map - (mean(map(~isnan(map))));     % Mean normalise
map(isnan(map))=0;                      % Set non-visited to zero 
ps=abs(fft2(map, nextpow2(size(map,1)), nextpow2(size(map,2))));      % zero padding the rate map



ps=fftshift(ps);                        % shifting low frequency components to the middle





cen=fftPad/2;
rCr=round((size(map,1)/2)*psCropFactor);  cCr=round((size(map,1)/2)*psCropFactor);
ps = ps( cen-rCr+1:cen+rCr, cen-cCr+1:cen+cCr );        % Crop power spectrum for display