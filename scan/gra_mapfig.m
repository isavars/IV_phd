
function [h_fig] = gra_mapfig(varargin)

% GRA_MAPFIG. Two methods of calling:
%
% 1) Call on command line
%
%           h_fig = gra_mapfig(data, map_ind, trials, cells);
%
% 2) Called by SCAN UI
%           
%           gra_mapfig();
%
%           - Called by base menu. 
% 			- Gets SCAN_DATA, evaluates .current_data, (will) look for default figure options.

% Notes on GRA_XXX function cascades:
% 			- Calls GRA_MULTIPLOT with correct proportions. This creates
% 			the figure window (but not the menus).
%             Store output (nice axis handle array) in figure 'userdata'.
% 			- Sets up custom menu bar. (Save as, Print, Options). (For callbacks, see below).
%             (Options Callback  - GRA_MAPFIG_UIOPTIONS(gcf);
% 			- Call GRA_MAPFIG_FILL. Pass handle to figure.
% 			- End


%%% Default Figure Options (This will come from SCAN_DATA default set, one day) %%%
opt.text_pos = 'tl';
opt.bm_thr = 0;
opt.n_steps = 11;
opt.scale_max = [];
opt.interp_uv = 1;
opt.handle = [];
opt.trials = [];
opt.trials_variable = 'default';
opt.cells = [];
opt.cells_variable = 'default';
opt.label = 'all';

%%% Check for calling mode, set relevant variables appropriately %%%
if ~isempty(varargin)
    % Called on command line %
    data = varargin{1};
    map_ind = varargin{2};
    opt.trials = varargin{3};
    opt.cells = varargin{4};
    fig_namer = data.rate_maps(map_ind).name;
    data_name = 'null';     % Can't always get original data name. Some options (set trial order) will not work.
    % Set options from input args %
    if length(varargin)>=6
        for ii=5:2:length(varargin)-1
            opt.(varargin{ii}) = varargin{ii+1};
        end
    end
        
else
    % Called by SCAN UI %
    SD = gss;
    data = evalin('base', SD.selData{1});
    % Check that a set of maps is selected %
    if isempty(SD.selMap)
        errordlg('Please select a set of maps.');
        return
    end
    opt.trials = SD.selTrial;
    opt.cells = SD.selCell;
    map_ind = SD.selMap;
    fig_namer = [SD.selData{1}, ' ', data.rate_maps(map_ind).name];
    data_name = SD.selData{1};
    
end

% Set up figure %
if min(size(data.rate_maps(map_ind).maps{1}))==1
    mapAR=1; % Map Aspect Ratio: set to one for directional maps ..
else
    % Otherwise set from first map %
    ii=1;
    while isempty(data.rate_maps(map_ind).maps{ii}) && ii<=numel(data.rate_maps(map_ind).maps)  % Need this protective loop in case first map is empty (from filtereing leaving no points, for examples)
        if ii==numel(data.rate_maps(map_ind).maps);  disp('All maps are empty!');   return;   end
        ii=ii+1;
    end
    mapAR=size(data.rate_maps(map_ind).maps{ii},2) / size(data.rate_maps(map_ind).maps{ii},1); 
end
h_fig = gra_multiplot(length(opt.cells), length(opt.trials), 'graphborder', repmat(0, [1 4]), 'figborder', [1 0.5 0.5 0.5],'plotsize',[2*mapAR 2]);
set(h_fig, 'defaulttextfontname', 'arial', 'numbertitle', 'off', ...
    'name', fig_namer);
set(h_fig,'inverthardcopy','off','color','white');      % Workaround for Matlab bug that prints white as text black with 'axis off'.

% Menus %
% set(h_fig, 'menubar', 'none');   % Turn off matlab default.
% uimenu('label', 'Save as eps ..', 'callback', 'scan_base(''saveEps'', gcf);');
uimenu('label', 'Print', 'callback', 'printdlg(gcf);');
uimenu('label', 'Options', 'callback', 'gra_mapfig_uioptions(gcf);');

%%% Assign figure properties (store in fig 'user') %%%
% Make dataset with maps + metadata only, for space saving %
temp_data.rate_maps = data.rate_maps(map_ind);
[temp_data.trials(1:length(data.trials)).user] = data.trials.user;
[temp_data.trials(1:length(data.trials)).trialname] = data.trials.trialname;
for ii=1:length(data.trials)
    [temp_data.trials(ii).cells(1:length(data.trials(ii).cells)).cellnum] = data.trials(ii).cells.cellnum;
    [temp_data.trials(ii).cells(1:length(data.trials(ii).cells)).tet] = data.trials(ii).cells.tet;
    [temp_data.trials(ii).cells(1:length(data.trials(ii).cells)).user] = data.trials(ii).cells.user;
end
% Assign to fig .userdata %
temp = get(h_fig, 'userdata');
temp.mapfig_options = opt;
%temp.map_ind = map_ind;
temp.data_name = data_name;
temp.data = temp_data;
set(h_fig, 'userdata', temp);

% Draw in maps %
gra_mapfig_fill(h_fig);