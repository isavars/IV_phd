function [fieldDistanceToWalls]=map_bordercell(map,ppm,binSize)
% Calculate border-score, B, as defined by Solstad et al, 2008.
%
%       B=map_bordercell(map)

% Important variables %
rateThresholdForFields = 0.5;   % Fraction of MR to use as field threshold.
sizeThresholdForFields = 200*0.36;  % This many cm-sq of contiguous bins.
e=5; % Edge trim factor. Ignore extreme-most e bins in visited env when defining edge.

% Find Fields %
hfPix=map > nanmean(map(:))*rateThresholdForFields;
fieldLabels = bwlabel(hfPix,4);
labelList=unique(fieldLabels);
fieldSizeThreshold = sizeThresholdForFields / ((100/ppm)*binSize);   % Number of bins to reach field size (200cmsq)
for ii=labelList'
    if sum(sum(fieldLabels==ii)) < fieldSizeThreshold
        fieldLabels(fieldLabels==ii) = 0;
    end
end

% Define wall positions %
[size_row, size_col] = size(map);
a = find(~isnan(map));           % find visited bins (for rows)
ainv = find(~isnan(map'));       % find visited bins (for columns)
if (e+1) > length(a);  e = length(a) - 1;  end   % To prevent crashing in case of very small visited enviroments.
wallPositions = [ceil(a(e+1)./size_col), ceil(a(end-e)./size_col), ...       % find edge columns
                 ceil(ainv(e+1)./size_row), ceil(ainv(end-e)./size_row)];    % find edge rows

% Calculate Cm, the maximum coverage of a field over the walls of the environment %





% Calculate Dm, the mean distance to the wall of all field pixels %
% First, find the minimum distance to a wall for each bin in the field %
[fieldRows fieldCols]=find(fieldLabels);
fieldCoords=[fieldRows fieldRows fieldCols fieldCols];          % 
wallPositionArray=repmat(wallPositions,size(fieldCoords,1),1);  % 
fieldDistanceToWalls=abs(fieldCoords - wallPositionArray);
minDistToWall=min(fieldDistanceToWalls,[],2);
% Next, weight these min distances to wall by firing rate %
fieldFRs=map(fieldLabels>0);       % The FR values should be indexed from the map in the same order as line "[fieldRows fieldCols]=find(fieldLabels);" above.
fieldFRs=fieldFRs./sum(fieldFRs);  % Normalise to total firing
minDistToWall=minDistToWall.*fieldFRs;






