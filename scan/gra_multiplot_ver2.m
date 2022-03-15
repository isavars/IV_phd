function [allHandles] = gra_multiplot_ver2(nRowTotal, nColTotal, varargin)
% Create a grid of multiple axes in a figure. If specified, rows will be split across more 
% than one figure window (so as to maintain a specified axes size, for example).
% 
%       handleStruct = gra_multiplot_ver2(R,C)
%       handleStruct = gra_multiplot_ver2(R,C, .. , options, .. );
% 
% Where R and C are the numbers of rows and columns of axes.
%
% Options:
%
% 'maxRowPerFig', N                     - Sets a maximum limit on how many rows in a figure window. Additional
%                                         axes after limit are created in new figure windows. Set to [] for no limit (default).
% 
% 'plotsize',[x,y]                      -  Sets the size of each plot. (Plot is graph axes plus border). Units 
%                                          are tenths of screen height (for both x and y). If figure is larger 
%                                          than the screen, plots are scaled down proportionally. Default = [2,2].
% 
% 'figborder', [left,right,top,bottom]  -  Sets the border around the edge of the figure. Units are tenths of 
%                                          screen height. Default = [0.5 0.5 0.5 0.5]. Borders are fixed even if graphs 
%                                          are scaled down.
%                                  
% 'axesborder',[left,right,top,bottom] -  Sets the border around the edge of each axis. Units are tenths of 
%                                          the plot size. Default = [1 1 1 1].
%
% OUTPUT:
%
%   handleStruct.axes       - Contains all axes handles, in a (R,C) array. The layout of axes in figure windows 
%                             is therefore not visible to any calling function.
%   handleStruct.rowLabels  - 1xR cell array, containing handles to row label text objects.
%   handleStruct.colLabels  - 1xC cell array, containing handles to col label text objects. For multiple figures, 
%                             each entry will be a vector, containing handles for text objects for specified column
%                             in each figure.
%   handleStruct.title      - 1xnFigure vector of handles to title text objects (titles are created in top-left corner
%                             of each figure).
%   handleStruct.figure     - 1xnFigure vector of handles to figures.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to make this function work for multiple figures (for the ultimate aim
% of having a 'max cells per figure' setting, or a 'lock plot size' option.
% 
% 1) gra_multiplot will become the function that handles the multiple figures. It will work out how many we need,
%    where they should be, and collect up all the important handles and pass them to the calling function, as if they
%    were all in one massive figure.
% 2) Existing gra_multiplot (to be renamed gra_subplot?) will exist pretty much as before, with two important changes:
%     2a) It will also create handles to the text labels as it goes, and pass the handles to these to gra_multiplot
%     2b) It will also have to accept an input for figure position
% 3) Note that the old model for passing handles (appdata property) will no longer work directly, as there is
%    more than figure. We can either put all handles under the first figure, or make the 

% 1) Plot size and max N plot  - Check how many figures needed. nPlotPerFig can be limited by either MaxNPlot, or screen_size, depending on whether
%                                maxNPlot*plotSize > screen_size (i.e. in a conflict, plotSize overrides maxNPlot). Then, nFig = totalPlot/nPlotPerFig.
%                                Call figures with locked plot size.
% 2) plot size only            - This is similar to 1, need to first calculate max N plot, on the basis of plot/screen.
% 3) Max N only                - Slightly different to 1 & 2. If N fig = 1 (when nPlot<maxNPLot), plot size will be free. Otherwise, will behave 
%                                as for 1&2. (need to calculate plot size)
% 4) No options                - Pass parameters directly to GRA_SUBPLOT


%%% Parse input %%%
% These are general defaults %
prms.figborder = ones(1,4) .* 0.05;
prms.axesborder = ones(1,4) .* 0.1;
prms.plotsize = [];                 % Leave this empty for now: will be set by user input or SCAn settings, or if not, will take value of prms.defaultPlotSize.
prms.defaultPlotSize = [0.2 0.2];   % This is the plotsize that will be used if there is no user specification.
prms.maxRowPerFig = [];           % Empty by default - means there is no maxPlot set, unless user specifified.
% If there is a fixed plot size or maxNPlot set in SCAn, this overrides the above defaults %
SD=gss;
if ~isempty(SD) && ~isempty(SD.settings.fixedPlotHeight)      % SD=[] if SCAn not running.
    prms.plotsize = ones(1,2).*SD.settings.fixedPlotHeight;   % SCAn setting refers to plot height only, and is one number.
    prms.maxRowPerFig = SD.settings.maxRowNumber;
end
% Any caller input to this function overrides these again %
for ii=1:2:length(varargin)
    prms.(varargin{ii}) = prms.(varargin{ii+1});
end

%%% Now set the plot size and number of rows per figure, depending what options are specified above %%%
if ~isempty(prms.plotsize) && ~isempty(prms.maxRowPerFig)
    % Case 1. Both plot size and max N plots specified. Both of these will be respected, although if plot_size*maxNplot > screen_size, the number of
    % plots per figure will be floor(screen_size/plot_size), not maxNplot.
    if prms.plotsize(2)*prms.maxRowPerFig > 1-sum(prms.figborder(3:4))   % Need to subtract figure border top and bottom.
        rowsPerFig = floor( 1-sum(prms.figborder(3:4)) / prms.plotsize(2));
    else
        rowsPerFig = prms.maxRowPerFig;
    end

elseif ~isempty(prms.plotsize) && isempty(prms.maxRowPerFig)
    % Case 2. Only plot size specified. Need to calculate rows per figure.
    rowsPerFig = floor( 1-sum(prms.figborder(3:4)) / prms.plotsize(2));
    
elseif isempty(prms.plotsize) && ~isempty(prms.maxRowPerFig)
    % Case (3). Max N only. Slightly different to 1 & 2. If N fig = 1 (when nPlot<maxNPLot), plot size will be global default. 
    % Otherwise, will behave as for 1&2. (need to calculate plot size). In all cases, assume square plot.
    if nRowTotal > prms.maxRowPerFig
        prms.plotsize = ones(1,2) .* (1-sum(prms.figborder(3:4)) / prms.maxRowPerFig);  % Assume square plot
        rowsPerFig = prms.maxRowPerFig;
    else
        prms.plotSize = prms.defaultPlotSize;
        rowsPerFig = nRowTotal;
    end
      
elseif isempty(prms.plotsize) && isempty(prms.maxRowPerFig)
    % Case 4. No parameters specified. Use global default plot size, and put all plots in one figure.
    rowsPerFig = nRowTotal;
    prms.plotSize = prms.defaultPlotSize;
    
end

% How many figures, and how many plots are on the last 'remainder' figure? %
nFigures = ceil(nRowTotal/rowsPerFig);
if nFigures>1
    remFigNRow = rem(nRowTotal,rowsPerFig);
    if remFigNRow>0
        nWholeFigures=nFigures-1;
    else
        nWholeFigures=nFigures;
    end
else
    nWholeFigures=1;
    remFigNRow = 0;
end

%%% Draw the figures, and collect the axes handles %%%
allHandles.axes = zeros(nRowTotal,nColTotal);
allHandles.rowLabels = cell(nRowTotal,1);   % colLabels need to be a cell array, so that each column label can be a vector of handles (representing the col titles
allHandles.colLabels = cell(nColTotal,1);   % on every figure), so that they can all be set using one function call. Package row labels the same, for consistency.
allHandles.titles = zeros(nWholeFigures+remFigNRow,1);
allHandles.figures = zeros(nWholeFigures+remFigNRow,1);
rowStartInd=1;
for ii=1:nWholeFigures
    hTemp=gra_subplot(rowsPerFig,nColTotal,'plotsize',prms.plotsize,'figborder',prms.figborder,'axesborder',prms.axesborder);
    rowInd=rowStartInd:rowStartInd+rowsPerFig;
    allHandles.axes(rowInd,:) = hTemp.axes;
    allHandles.rowLabels{rowInd} = deal(hTemp.rowLabels);
    for jj=1:nColTotal
        allHandles.colLabels{jj} = [allHandles.colLabels{jj} hTemp.colLabels(jj)];
    end
    allHandles.title(ii) = hTemp.title;
    allHandles.figure(ii) = hTemp.figure;
    rowStartInd=rowStartInd + rowsPerFig;
end
% .. this is for the 'remainder' figure, if there is one %
if remFigNRow>0
    hTemp=gra_subplot(remFigNRow,nColTotal,'plotsize',prms.plotsize,'figborder',prms.figborder,'axesborder',prms.axesborder);
    rowInd=rowStartInd:rowStartInd+remFigNRow;
    allHandles.axes(rowInd,:) = hTemp.axes;
    allHandles.rowLabels{rowInd} = deal(hTemp.rowLabels);
    for jj=1:nColTotal
        allHandles.colLabels{jj} = [allHandles.colLabels{jj} hTemp.colLabels(jj)];
    end
    allHandles.title(end) = hTemp.title;
    allHandles.figure(end) = hTemp.figure;
end
    


