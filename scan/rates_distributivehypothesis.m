function [dirDHMap placeDHMap]=rates_distributivehypothesis(pos,spk,dirSm,placeSm)
% Create distributive hypothesis directional heading map.
%   [dirDHMap placeDHMap]=rates_distributivehypothesis(pos,spk,dirSmooth,placeSmooth)
% See Muller, Bostock, Taube & Kubie, 1994.

warning('off', 'MATLAB:divideByZero');
[nX nY nD]=size(pos);

% Make rate maps %
placeRateMap=sum(spk,3)./sum(pos,3);
placeRateMap=rates_smooth(placeRateMap,placeSm,isnan(placeRateMap));
dirRateMap=squeeze(sum(sum(spk,1),2)) ./ squeeze(sum(sum(pos,1),2));
dirRateMap=rates_smooth(dirRateMap,dirSm);

warning('off', 'MATLAB:divideByZero'); % RATES_SMOOTH sets warning state to on. Neat to fix this.

% Generate place x dir dwell time array (2D, place index is linearised) %
pos=permute(pos,[3,1,2]);
pxdDwell=reshape(pos,[nD nY*nX]);   % This is in format row=dir bin, column=all place bins

% Expected rates: direction %
placeRates=reshape(placeRateMap,[nY*nX 1]);                      % This is the rate map, linearised as a column vector, same order as place bins in pxdDwell
expSpkDir = nansum(  pxdDwell' .* repmat(placeRates,[1 nD])  );  % Expected number of spikes. Sum, for each dir bin, the dir dwell times over all place bins, weighted by firing rate in each place bin.
dirPosMap=sum(pxdDwell,2);
dirDHMap = ( expSpkDir' ./ dirPosMap );         % Dist hypo rate map: expected number of spikes, / by dir dwell times.

% Expected rates: place %
expSpkPl = nansum(  pxdDwell .* repmat(dirRateMap,[1 nX*nY])  );  % Expected number of spikes. Sum, for each place bin, the place dwell times over all dir bins, weighted by firing rate in each dir bin.
plPosMap = sum(pxdDwell,1);
placeDHMap = ( expSpkPl./plPosMap )';
placeDHMap = reshape(placeDHMap,[nY nX]);

% Smooth final output %
dirDHMap = rates_smooth(dirDHMap,dirSm);
placeDHMap = rates_smooth(placeDHMap,placeSm,isnan(placeRateMap));

warning('on', 'MATLAB:divideByZero');






    





