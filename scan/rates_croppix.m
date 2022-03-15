function [xCrop,yCrop] = rates_croppix(x,y,coords,e)
% Find centre and crop
% 
%       [x,y] = rates_crop(x,y,coords,e);
%
% coords format:
%  [A B]     - Crop a box, AxB (width x height), around the environment centre.
%  [x y A B] - Crop absolute co-ordinates, defined [top_left_x,top_left_y,width,height].
%
% e = Edge trimming variable. Set 0 to ?. Not relevant for absolute coords,
% but need to give [] as placeholder.

if length(coords)==4
    %%% Crop in absolute space - just subtract top-left corner from x and y %%%
    xCrop=x-coords(1);
    yCrop=y-coords(2);
elseif length(coords)==2
    %%% Crop centred on environment %%%
    % Find edges %
    xSort=sort(data.trials(jj).x);
    xExtent=[xSort(e) xSort(end-e)];
    ySort=sort(data.trials(jj).y);
    yExtent=[ySort(e) ySort(end-e)];
    % Find centre and edges %
    cen=floor([mean(xExtent), mean(yExtent)]);
    xCropEdge=cen(1) - round(coords(1)/2);
    yCropEdge=cen(2) - round(coords(2)/2);
    xCrop=x-xCropEdge;
    yCrop=y-yCropEdge;
end