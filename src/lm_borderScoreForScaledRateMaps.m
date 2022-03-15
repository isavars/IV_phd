function [borderScore]=lm_borderScoreForScaledRateMaps(map,binSize,varargin)
% Calculate border-score, B, as defined by Solstad et al, 2008.
%
%       B=map_bordercell(map,binSizeInCmSq);
%       B=map_bordercell(map,binSizeInCmSq, .. options, .. );
%
% Options:
%
%   'rateThr',  0.2,    - Threshold for defining fields, proportion of peak rate.
%   'sizeThr',  200,    - Threshold for minimum field size to enter into the analysis.


% Important variables %
prms.rateThr = 0.3;       % Fraction of max rate to use as field threshold.
prms.sizeThr = 200;       % This many cm-sq of contiguous bins.
for ii=1:2:length(varargin)
    prms.(varargin{ii}) = varargin{ii+1};
end


% Find Fields %
bwMap = map >= ( nanmax(map(:))*prms.rateThr );
% bwMap = map >= 1;
fieldLabels = bwlabel(bwMap,4);
labelList = setdiff(unique(fieldLabels),0);
% Cycle through fields, if a field is big enough, get its wall coverage (hence to calculate Cm) %
fieldsForDm = zeros(size(map));
coverageOnFourWalls = nan(length(labelList),4);
for ii=1:length(labelList)
	fieldTemp = fieldLabels==labelList(ii);
    if (sum(fieldTemp(:)) * binSize) < prms.sizeThr   % If field too small, continue
        continue
    end
    % Get coverage of single fields over individual walls %
    coverageOnFourWalls(ii,1) = sum(fieldTemp(:,1));
    coverageOnFourWalls(ii,2) = sum(fieldTemp(:,end));
    coverageOnFourWalls(ii,3) = sum(fieldTemp(1,:));
    coverageOnFourWalls(ii,4) = sum(fieldTemp(end,:));
    % Keep a record of the position of this field, so we have all 'big enough' fields for calculation of Dm, below %
    fieldsForDm = fieldsForDm | fieldTemp;
end
% If no big enough fields, return %
if sum(fieldsForDm(:)) == 0
    borderScore = nan;
    return
end
% Calculate Cm, the maximum coverage of a field over the walls of the environment %
Cm = nanmax(coverageOnFourWalls(:)) / size(map,1);

% Calculate Dm, the mean distance to the wall of all field pixels, weighted by firing rate %
% Make a matrix showing distance to nearest wall %
hMapSize = ceil(size(map,1)/2);
if rem(size(map,1),2)
    [Y,X]=meshgrid([1:hMapSize fliplr(1:hMapSize-1)]);
else
    [Y,X]=meshgrid([1:hMapSize fliplr(1:hMapSize)]);
end
dist2Wall = squeeze(min(  cat(3,X,Y)  ,[],  3  ));
% Multiply this by the normalised FR map, and then get the mean distance to wall for the 'big enough' field pixels %
rateInFields = map(fieldsForDm);
dist2WallInFields = dist2Wall(fieldsForDm);
Dm = transpose(dist2WallInFields) * (rateInFields / sum(rateInFields));
Dm = Dm / (size(map,1)/2);

% Calculate the border score, (Cm-Dm) / (Cm + Dm) %
borderScore = (Cm-Dm) / (Cm + Dm);







