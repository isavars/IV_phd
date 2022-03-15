function [rtn]=map_popvectcrosscorr(mapsA,mapsB,varargin)
% Population vector-based spatial cross-correlogram.
%
%       C=map_popvectcrosscorr(mapsA,mapsB)
%
% Functions similarly to MAP_POPVECT, but also offsets mapsB by varying
% amounts, so as to create a spatial cross-correlogram.
%
% Can also rotate one set of maps to search for the best fit:
%
%       C=map_popvectcrosscorr(mapsA,mapsB,rotation_angles)
%
% See MAP_POPVECT and Fyhn, Hafting et al, 2007, Nature, for more info.

if isempty(varargin)
    rotAngles = 0;
else
    rotAngles = varargin{1};
end

stackA = cat(3,mapsA{:});
stackB = cat(3,mapsB{:});

for aa=1:length(rotAngles)
    % Rotate map B (will perform one rotation of 0deg by default %
    visMask = ~isnan(stackB);                              %
    stackBRot=imrotate(stackB,rotAngles(aa),'crop');       % Note that, when rotating, need to also keep track of 
    visMaskRot=imrotate(visMask,rotAngles(aa),'crop');     % unvisited bins and reset to NaN.
    stackBRot(~visMaskRot) = NaN;                          %

    %%% This is the beginning of the core Offset + Pop Vect Function %%%
    nRows=size(stackA,1);
    nCols=size(stackA,2);
    offsetsRows=1:(nRows*2);    % These are the lists of row and column offsets.
    offsetsCols=1:(nCols*2);    % (2 times length of original map in both directions)
    crossCorr = nan( nRows*2, nCols*2 );

    for ii=offsetsRows
        parfor jj=offsetsCols

            % Create each instance of the offset map %
            offsetStack = nan(size(stackA,1)*3,size(stackA,2)*3,size(stackA,3));  % Pre-allocate offset matrix, 3x3 grid of original (but all NaN).
            offsetStack(ii:ii+nRows-1, jj:jj+nCols-1, :) = stackBRot;             % Copy in stackB at given offset from top-left corner of 3x3.
            offsetStack = offsetStack(nRows+1:nRows*2, nCols+1:nCols*2, :);       % Cut out 'centre' stack area (middle of 3x3 grid), for final offset maps.

            % Create a mask of mutually visited bins, check this is large enough to run pop vect %
            visMaskA = sum(isnan(stackA),3) == 0;
            visMaskB = sum(isnan(offsetStack),3) == 0;
            visMask = visMaskA & visMaskB;
            if sum(visMask(:)) < 20
                continue
            end
            [visRow visCol] = find(visMask);

            % Run Pop Vect Analysis %
            dotProd=nan(1,sum(visMask(:)));
            for kk=1:length(dotProd)
                a = double(squeeze( stackA(visRow(kk),visCol(kk),:) ));
                b = double(squeeze( offsetStack(visRow(kk),visCol(kk),:) ));
                dotProd(kk) = dot(a,b) / (norm(a)*norm(b));
            end
            crossCorr(ii,jj)=mean(dotProd); 

        end
    end
    rotCrossCorrs{aa}=crossCorr; %#ok<AGROW> % This is the final crossCorr for rotation aa.
end
    
if length(rotAngles)==1
    rtn=rotCrossCorrs{1};
else
    rtn=rotCrossCorrs;
end
    
    
    
%%%%%%%%%%%%%% Scaling code - removed %%%%%%%%%%%%%%%%%%%%%%%%
% scaleFactors = 0.8:0.1:1.2;
% hFig=gra_multiplot(2,length(scaleFactors)); hAxArray = getappdata(hFig,'axesHandles');
% mapsBPreserved = mapsB;
%     
% for aa=1:length(scaleFactors)
% 
%     mapsB = rates_transform(mapsBPreserved,1,mapsA{1},1,'force_scale',scaleFactors(aa));
%     gra_plotmap(mapsB{end},'handle',hAxArray(1,aa));
%     gra_plotmap(mapsA{end},'handle',hAxArray(2,aa));
%     drawnow
%     
%     scaleCrossCorrs{aa} = crossCorr;
% end
%     
% rtn=scaleCrossCorrs;
    