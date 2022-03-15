
function [H] = gra_multilabel(fh, mode, labels, varargin)

% Label the rows or columns of a multiPlot created figure.
%
%   H = gra_multilabel(figure_handle, mode, labels);
%   H = gra_multilabel(figure_handle, mode, labels, .. 'options' .. );
%
% figure_handle is the handle to a gra_multiplot plot, mode is 'row' or 'col' to label row or column respectively,
% labels is a cell array of strings containing the labels. 
%
% If mode is 'title', will place a title in the top left-hand corner. Labels argument should
% be a string, or a cell array of strings for multi-line titles.
%
% Options are any valid Matlab text property name and value pairs, given in the standard format.
%
% H is an array of handles to the resulting text objects

temp = get(fh, 'userdata'); 
h_array = temp.axes_handle_array;

% Sort out the format of 'labels'. %
if iscategorical(labels)
    for ii=1:length(labels);   new{ii} = char(labels(ii));    end
    labels = new;
elseif isnumeric(labels)
    for ii=1:length(labels);   new{ii} = num2str(labels(ii));    end
    labels = new;
elseif islogical(labels)
    for ii=1:length(labels);   new{ii} = num2str(double(labels(ii)));    end
    labels = new;
end

%%% Draw Labels %%%
switch mode
case 'row'
    % Label Rows %
    n_plots = size(h_array, 1);
    
    for i=1:n_plots
        if i>length(labels)
            break
        end

        axes(h_array(i,1));
        
        axes_pos = get(gca, 'position');
        border_width = axes_pos(1) / axes_pos(3);
        
        if size(h_array, 2) > 1                              % Correct for graph border. Needs to look at next column to
            axes2_pos = get(h_array(i,2), 'position');       % the right, so will only do if this exists.
            graph_border_corr = (axes2_pos(1) - axes_pos(1) - axes_pos(3)) / 2;
            border_width = border_width + graph_border_corr;
        end
        
        text_pos = (border_width / 2) * -1;
        
        if isnumeric(labels{i})
            labels{i} = num2str(labels{i});
        end
        H(i) = text('units', 'normalized', 'position', [text_pos, 0.5], 'string', labels{i}, 'horizontalalignment', 'center', varargin{1:end});
        
    end
case 'col'
    % Label columns %
    n_plots = size(h_array, 2);

    for i=1:n_plots
        if i>length(labels)
            break
        end
        
        axes(h_array(1,i));
        
        axes_pos = get(gca, 'position');
        border_width = (1 - (axes_pos(2) + axes_pos(4))) / axes_pos(4);
        text_pos = (border_width / 2) + 1;
        
        if isnumeric(labels{i})
            labels{i} = num2str(labels{i});
        end
        H(i) = text('units', 'normalized', 'position', [0.5, text_pos], 'string', labels{i}, 'horizontalalignment', 'center', varargin{1:end});
        
    end
case 'title'
        % Put a title in the top left-hand corner %
        axes(h_array(1,1));
        axes_pos = get(gca, 'position');
        
        % position on left %
        border_width = axes_pos(1) / axes_pos(3);
        if size(h_array, 2) > 1                              % Correct for graph border. Needs to look at next column to
            axes2_pos = get(h_array(1,2), 'position');       % the right, so will only do if this exists.
            graph_border_corr = (axes2_pos(1) - axes_pos(1) - axes_pos(3)) / 2;
            border_width = border_width + graph_border_corr;
        end
        text_pos_left = (border_width * 0.9) * -1;
        
        % Position at the top %
        border_width_top = (1 - (axes_pos(2) + axes_pos(4))) / axes_pos(4);
        text_pos_top = (border_width_top / 2) + 1;
        
        H = text('units', 'normalized', 'position', [text_pos_left, text_pos_top], 'string', labels, 'horizontalalignment', 'left', 'fontweight', 'bold', varargin{1:end});
end
    
