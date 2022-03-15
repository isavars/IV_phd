function [data,pathLimCheck] = lm_meanRateMap_findEdgesAndScalePath_CA1(data,minOccForEdge,boxExtent)
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

prms.binSize = 10;

pathLimCheck = nan(length(data.trials),6);
for ii=1:length(data.trials)
    
    if strcmp(data.trials(ii).user.environment,'hp_large')
        boxExtentTrial = boxExtent*2;
        %%LM edit: there is some adult data set that are recorded in circle as novel_light
        %I tagged these in SCAN - circle was diam:80cm, so 80/2.5 = 32 bins
    elseif strcmp(data.trials(ii).user.environment,'novel_light') && isfield(data.trials(ii).user,'typeNovel')
        boxExtentTrial = 320;
    else
        boxExtentTrial = boxExtent;
    end
    
    dims = {'x','y'};
    for jj=1:2
        
        path = double(  data.trials(ii).( dims{jj} )  );
        
        %%LM: this is to cope with data sets which contain reflections -
        %%have to deal with individually; 18/04/15: have checked that this
        %%is the only trial/dataset where reflections cause problem for path
        %%scaling
        if isfield(data,'load_selection') && strcmp(data.load_selection.dataName,'R1783_100924') && strcmp(data.trials(ii).user.environment,'novel_dark')
            
            path_binned = floor(path./prms.binSize);
            binMin = min(path_binned);
            maxPoss = binMin+round(boxExtent/prms.binSize);
            ind2remove = path_binned > maxPoss;
            path(ind2remove) = NaN;
     
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        pathHist = histc( path, 0.5:1:768.5  );    % 768 is the maximum extent for the DACQ camera.
        
        lowerEdge = find(pathHist>=minOccForEdge, 1, 'first');
        upperEdge = find(pathHist>=minOccForEdge, 1, 'last');
        
        path( path>upperEdge ) = NaN;
        path( path<=lowerEdge ) = NaN;       % Doing <=lowerEdge, then subtracting lowerEdge (line 23), makes the lower limit zero, and therefore the first pixel 1.something.
        
        path = path-lowerEdge;
        path = path.*(  boxExtentTrial / (upperEdge-lowerEdge) );
        
        path(path>boxExtentTrial) = boxExtentTrial;    % Have checked that when this happens, it is a many decimal places rounding error.
        
        pathLimCheck( ii, ((jj-1)*3)+1  ) = min(path);
        pathLimCheck( ii, ((jj-1)*3)+2  ) = max(path);
        pathLimCheck( ii, ((jj-1)*3)+3  ) = upperEdge - lowerEdge;
        
        data.trials(ii).( dims{jj} ) = path;
        
    end
    % For consistency, make sure that all x=nan and all y= nan match up %
    data.trials(ii).x(  isnan(data.trials(ii).y)  ) = nan;
    data.trials(ii).y(  isnan(data.trials(ii).x)  ) = nan;
    
    
end