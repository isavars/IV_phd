
function [] = gra_mapfig_fill(h_fig)

% Draw all maps into figure.
%
%       gra_mapfig_fill(h_fig);
%
% Original notes:
%
% GRA_MAPFIG_FILL	- Mainly just a loop to draw all maps into axes, calling GRA_DRAWMAP with correct options.
% 			        - Needs to know 'order' and 'ori' to know where to put maps.
% 			        - Could also be responsible for labelling axes.

%%% NOTE, 18/12/06: Map drawing with multiple selections - the temporary workaround is to draw the first map set only (see line 30)

temp = get(h_fig, 'userdata');
data = temp.data;
opt = temp.mapfig_options;
h_array = temp.axes_handle_array;
series_ind = 1;

%%% Draw maps into axes, rearranging axes handles array %%%
for ii = 1:length(opt.trials)        
    for jj = 1:length(opt.cells)
        % Draw maps into selected axis %
        axes(h_array(jj, ii));
        gra_plotmap( double(  data.rate_maps(1).maps{ opt.trials(ii),opt.cells(jj),series_ind }  ), opt);     
    end
end

%%% Label %%%
if strcmp(opt.label, 'trial') || strcmp(opt.label, 'all')
    gra_multilabel(h_fig, 'col', getdatanames(data, 'trial', 'short', opt.trials));
end
%%% TODO: if you want labels that include user tags, probably best way to do it is to loop through data here, and make nested
%%%       cell array of '{trialName','tag'} for each trial.
if strcmp(opt.label, 'cell') || strcmp(opt.label, 'all')
    gra_multilabel(h_fig, 'row', getdatanames(data, 'cell', 'short', opt.cells));
end
if ~strcmp(temp.data_name,'null')
    gra_multilabel(h_fig,'title',temp.data_name);
end
