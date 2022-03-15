
function [h_fig] = gra_multiplot(n_r, n_c, varargin)

% Draw a figure containing multiple axes.
% 
%       fig_handle = gra_multiplot(R,C)
%       fig_handle = gra_multiplot(R,C, .. , options, .. );
% 
% Where R and C are the numbers of rows and columns of graphs. 
%
% To get the handles of the axes created, run:
%
%   A = getappdata(fig_handle,'axesHandles');
%
% The handle in A(r,c) matches the axis at row=r, col=c in the figure.
%
% Options:
% 
% 'plotsize',[x,y]                      -  Sets the size of each plot. (Plot is graph axes plus border). Units 
%                                          are tenths of screen height (for both x and y). If figure is larger 
%                                          than the screen, plots are scaled down proportionally. Default = [2,2].
% 
% 'figborder', [left,right,top,bottom]  -  Sets the border around the edge of the figure. Units are tenths of 
%                                          screen height. Default = [0.5 0.5 0.5 0.5]. Borders are fixed even if graphs 
%                                          are scaled down.
%                                  
% 'axesborder',[left,right,top,bottom]  -  Sets the border around the edge of each axis. Units are tenths of 
%                                          the plot size. Default = [1 1 1 1].
%
% 'axOrPanel', {'axes', 'panel'}        -  After the figure is split into the grid, do we insert an axis or a uipanel object at each grid point?
%                                          Default is 'axes'.
% 'polAx',     [rowIndices; colIndices] -  LM edit; added possibilty of defining polar axes (for use of 'polarhistogram'); just supply all indices of axes that are supposed to be pax object

                                 

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

objectInGrid = 'axes';
polAxInd = false(n_r, n_c);

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
                case 'axesOrPanel'
                    objectInGrid = varargin{i+1};
                case 'polAx'
                    ind = sub2ind([n_r n_c], varargin{i+1}(1,:), varargin{i+1}(2,:));
                    polAxInd(ind) = true;
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
if fig_width*(1/screen_proportion) > fig_height
    paper = 'landscape';
else
    paper = 'portrait'; 
end

    
corner_x = (1 - fig_width) / 2;     % Set bottom left corner of figure
corner_y = (1 - fig_height) / 2;    %

%%% Draw figure %%
h_fig = figure('units', 'normalized', 'position', [corner_x, corner_y, fig_width, fig_height], ...
               'paperorientation', paper, 'paperpositionmode', 'auto','defaulttextinterpreter','none','toolbar','none',...
               'inverthardcopy','off','color','white','papertype','A4','renderer','painters');
uimenu('label', 'Save as eps ..', 'callback', 'scan_base(''saveEps'', gcf);');
uimenu('label', 'Save as pdf ..', 'callback', 'scan_base(''savePdf'', gcf);');


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
        
        if strcmp( objectInGrid, 'axes' ) && ~polAxInd(j,i)
            axes_array(j,i) = axes('position', [((i-1)*plot_width) + graph_margin_width + figure_margin_width, (1-(j*plot_height)) + graph_margin_height - figure_margin_height, ...
                axes_width, axes_height], 'xticklabel', [], 'yticklabel', [], 'defaulttextinterpreter', 'none');
        elseif polAxInd(j,i)
            axes_array(j,i) = polaraxes('position', [((i-1)*plot_width) + graph_margin_width + figure_margin_width, (1-(j*plot_height)) + graph_margin_height - figure_margin_height, ...
                axes_width, axes_height], 'rticklabel', [], 'thetaticklabel', [], 'ticklabelinterpreter', 'none');
        else
            axes_array(j,i) = uipanel('position', [((i-1)*plot_width) + graph_margin_width + figure_margin_width, (1-(j*plot_height)) + graph_margin_height - figure_margin_height, ...
                axes_width, axes_height]);
        end
        
    end
end

%%% Put nice layout axis handles in user data for figure. %%
temp.axes_handle_array = axes_array;
set(gcf, 'userdata', temp);
setappdata(h_fig,'axesHandles',axes_array);




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