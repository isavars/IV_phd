function [rtn] = gra_subplot(n_r, n_c, varargin)
% Create a figure containing a grid of multiple axes.
% 
%       handleStruct = gra_subplot(R,C)
%       handleStruct = gra_subplot(R,C, .. , options, .. );
% 
% Where R and C are the numbers of rows and columns of graphs. 
%
% OUTPUT HANDLE STRUCT:
%
%     handleStruct.axes       - Handles to axes, in an (R,C) array.
%     handleStruct.rowLabels  - Handles to row label text objects (one for each row).
%     handleStruct.colLabels  - Handles to column label text objects (one for each column).
%     handleStruct.title      - Handle to title text object (in top-left of figure border).
%     handleStruct.figure     - Figure Handle
%
% INPUT OPTIONS:
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



%%% DEFAULT %% shows where there is a default setting for figure appearance, that can be altered a user input option.

plot_height_screen = 0.2;   %% DEFAULT %% this sets default height of a plot, in relation to screen size ..
plot_width_screen = 0.2;    %% DEFAULT %%  .. units are normalized screen height units (whole screen = 1), for both the width and height of plot.

figure_border_left = 0.05;       %
figure_border_right = 0.05;      %
figure_border_bottom = 0.05;     %
figure_border_top = 0.05;        %

graph_border_left = 0.1;        % 
graph_border_right = 0.1;       %
graph_border_top = 0.1;         % 
graph_border_bottom = 0.1;      %


%%% Check user input for options %%
if ~isempty(varargin)
    for i=1:length(varargin)
        if ischar(varargin{i})
            % Always / 10, because user units are tenths of screen, tenths of plot %
            switch varargin{i}
            case 'plotsize'
                plot_height_screen = (varargin{i+1}(2)) / 10;
                plot_width_screen = (varargin{i+1}(1)) / 10;
            case 'figborder'
                figure_border_left = (varargin{i+1}(1)) / 10;
                figure_border_right = (varargin{i+1}(2)) / 10;
                figure_border_top = (varargin{i+1}(3)) / 10;
                figure_border_bottom = (varargin{i+1}(4)) /10;
%             case 'figborderleft'
%                 figure_border_left = varargin{i+1} / 10;
%             case 'figborderright'
%                 figure_border_right = varargin{i+1} / 10;
%             case 'figbordertop'
%                 figure_border_top = varargin{i+1} / 10;
%             case 'figborderbottom'
%                 figure_border_bottom = varargin{i+1} / 10;
            case 'graphborder'  % This nomenclature is left functioning so as not to break old functions.
                graph_border_left = (varargin{i+1}(1)) / 10;
                graph_border_right = (varargin{i+1}(2)) / 10;
                graph_border_top = (varargin{i+1}(3)) / 10;
                graph_border_bottom = (varargin{i+1}(4)) / 10;
            case 'axesborder'  % 'axesborder' is the new 'graphborder', for consistency with matlab terminology.
                graph_border_left = (varargin{i+1}(1)) / 10;    
                graph_border_right = (varargin{i+1}(2)) / 10;
                graph_border_top = (varargin{i+1}(3)) / 10;
                graph_border_bottom = (varargin{i+1}(4)) / 10; 
%             case 'graphborderleft'
%                 graph_border_left = varargin{i+1} / 10;
%             case 'graphborderright'
%                 graph_border_right = varargin{i+1} / 10;
%             case 'graphbordertop'
%                 graph_border_top = varargin{i+1} / 10;
%             case 'graphborderbottom'
%                 graph_border_bottom = varargin{i+1} / 10;
%             otherwise
%                 error(['multiPlot: layout parameter ''', varargin{i}, ''' not recognized.']);
            end
        end
    end
end


%%% Correct for screen rectangle, so images with specified as squares display as squares %%
screen_pix = get(0, 'screensize');
screen_proportion = screen_pix(4) / screen_pix(3);

plot_width_screen = plot_width_screen * screen_proportion;          % 
figure_border_right = figure_border_right * screen_proportion;      % Correct all horizontal measurements
figure_border_left = figure_border_left * screen_proportion;        %


%%% Set figure size. If figure goes off screen, reduce proportionally %%%
figure_border_width = figure_border_left + figure_border_right;
figure_border_height = figure_border_top + figure_border_bottom;
figure_plots_width = plot_width_screen * n_c;
figure_plots_height = plot_height_screen * n_r;

% Because figure borders are set in terms of screen proportion units, size limit test looks at total plot sizes %
% within borders, then reduces these. Total figure size is calculated later, below, by adding borders to plots  %
if figure_plots_width > (0.9 - figure_border_width)
    figure_plots_height = figure_plots_height * ((0.9 - figure_border_width) / figure_plots_width);
    figure_plots_width = 0.9 - figure_border_width;  
end
if figure_plots_height > (0.9 - figure_border_height)
    figure_plots_width = figure_plots_width * ((0.9 - figure_border_height) / figure_plots_height);
    figure_plots_height = 0.9 - figure_border_height;  
end

fig_height = figure_plots_height + figure_border_height;
fig_width = figure_plots_width + figure_border_width;

% Get paper mode here %
if fig_width > fig_height; paper = 'landscape'; else paper = 'portrait'; end;
    
corner_x = (1 - fig_width) / 2;     % Set bottom left corner of figure
corner_y = (1 - fig_height) / 2;    %

%%% Draw figure %%
h_fig = figure('units', 'normalized', 'position', [corner_x, corner_y, fig_width, fig_height], ...
               'paperorientation', paper, 'paperpositionmode', 'auto','defaulttextinterpreter','none','toolbar','none',...
               'inverthardcopy','off','color','white','papertype','A4');
uimenu('label', 'Save as eps ..', 'callback', 'scan_base(''saveEps'', gcf);');


%%% Set plot + figure margin size (Units are normalized figure) %%%
border_proportion_width = figure_border_width / fig_width;
border_proportion_height = figure_border_height / fig_height;

plot_width = (1 - border_proportion_width) / n_c ;      % This is the (hypothetical) space given to
plot_height = (1 - border_proportion_height) / n_r;     % each graph, including borders.

figure_margin_width = figure_border_left / fig_width;          % Sets the margin of the figure
figure_margin_height = figure_border_top / fig_height;         % around the left/top of the figure

%%% Set graph + graph margin size (Units are normalized plot) %%%          
axes_width = plot_width * (1 - graph_border_left - graph_border_right);
axes_height = plot_height * (1 - graph_border_top - graph_border_bottom);
graph_margin_width = plot_width * graph_border_left;                                % Sets the position of the graph axes in the plot.
graph_margin_height = plot_height * graph_border_bottom;


%%% Create Axes %%%
for i=1:n_c     
	for j=1:n_r
        axes_position = [((i-1)*plot_width) + graph_margin_width + figure_margin_width, (1-(j*plot_height)) + graph_margin_height - figure_margin_height, axes_width, axes_height];
        axes_array(j,i) = axes('position',axes_position, 'xticklabel', [], 'yticklabel', [], 'defaulttextinterpreter', 'none'); 
    end
end

%%% Create text objects for labels, in the middle of the left and top figure borders %%%
% (Do this by selecting the left- and top-most lines of axes, and creating text objects that are children of those axes,   %
% in the space outsode the axis where the figure border is. 'title' is a single text at the top-left corner of the figure. %
% Row labels %
leftMarginInAxesCoords=(figure_margin_width/axes_width) + graph_margin_width;
xPos= -leftMarginInAxesCoords/2;
for ii=1:n_r
	hRow(ii)=text('parent',axes_array(ii,1),'units','normalized','position',[xPos, 0.5],'string','','horizontalalignment','center');
end
% Column labels %
topMarginInAxesCoords=(figure_margin_height/axes_height) + graph_margin_height;
yPos= -topMarginInAxesCoords/2;
for ii=1:n_c
	hCol(ii)=text('parent',axes_array(1,ii),'units','normalized','position',[0.5,yPos],'string','','horizontalalignment','center');
end
hTitle=text('parent',axes_array(1,1),'units','normalized','position',[xPos,yPos],'string','','horizontalalignment','center''fontweight','bold');


%%% Put nice layout axis handles in return struct. %%
rtn.axes=axes_array;
rtn.figure=h_fig;
rtn.rowLabels=hRow;
rtn.colLabels=hCol;
rtn.title=hTitle;







%%%%% Junk %%%%%%%%

% This was deleted from help. Functionality is still there, just I was never using it
% and it was cluttering up the help section.
%
% 'figborderleft', value            -   Sets each one of the figure borders individually
% '    "    right', value
% '    "    top', value
% '    "    bottom', value
%
%'graphborderleft',           -   Set each graph border individually.
% '     "     right',
% '     "     top',
% '     "     bottom',