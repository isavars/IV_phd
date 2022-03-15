function [data,pathLimCheck] = lm_meanRateMap_findEdgesAndScalePath_v2(data,minOccForEdge,boxExtent,envCheck)
% Find the edges of vis env, and scale path inside reamining vis env so that it is exactly
% equivalent to 62.5cm in each dimension (so that each binned env is exactly 250 pixels across at 
% 400ppm, therefore 25bins when using 2.5cm bins).
%   
%       [data,pathLimCheck] = lm_meanRateMap_findEdgesAndScalePath(data,minOccForEdge,boxExtent)
%
% minOccForEdge sets the threshold number of pos samples for the edge of the visted environment.
% boxExtent is the size (side length) of the box after scaling.
%
% pathLimCheck shows (by row):  [scaled lower x, scaled upper x, non-scaled x extent, scaled lower y, scaled upper y, non-scaled y extent]. One
% row shows numbers from one trial.


pathLimCheck = nan(length(data.trials),6);
envCheck = {'fam','nov','sleep'}; %IV 11/02/2020
for ii=1:length(data.trials)
    
    % skip if not correct environment
    if ~any( strcmp(data.trials(ii).user.environment,envCheck) )  % TW edit: ~any(   ) allows envCheck to also be a cell array list of envs to scale.
        continue
    end
    
    dims = {'x','y'};
    for jj=1:2
        
        path = double(  data.trials(ii).( dims{jj} )  );
        
        pathHist = histc( path, 0.5:1:768.5  );    % 768 is the maximum extent for the DACQ camera.
        
        lowerEdge = find(pathHist>=minOccForEdge, 1, 'first');
        upperEdge = find(pathHist>=minOccForEdge, 1, 'last');
        
        path( path>upperEdge ) = NaN;
        path( path<=lowerEdge ) = NaN;       % Doing <=lowerEdge, then subtracting lowerEdge (line 23), makes the lower limit zero, and therefore the first pixel 1.something.
        
        path = path-lowerEdge;
        path = path.*(  boxExtent / (upperEdge-lowerEdge) );
        
        path(path>boxExtent) = boxExtent;    % Have checked that when this happens, it is a many decimal places rounding error.
        
        pathLimCheck( ii, ((jj-1)*3)+1  ) = min(path);
        pathLimCheck( ii, ((jj-1)*3)+2  ) = max(path);
        pathLimCheck( ii, ((jj-1)*3)+3  ) = upperEdge - lowerEdge;
        
        data.trials(ii).( dims{jj} ) = path;
        
    end
    % For consistency, make sure that all x=nan and all y= nan match up %
    data.trials(ii).x(  isnan(data.trials(ii).y)  ) = nan;
    data.trials(ii).y(  isnan(data.trials(ii).x)  ) = nan;
    
    
end