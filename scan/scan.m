
function [] = scan(varargin)

% Open SCAn (Spatial Correlate Analysis)
%
% Just type 'scan'

% This is the SCAn starting function:
% - checks if there is an instance already running, and if not
% - opens GUI base fig file
% - initialises app data struct and saves to gui.
% - Puts user-defined menus on base window (calls SCAN_BASE_MENUS)
% - calls SCAN_BASE('listRefreshData', hListData) , to populate window.

% Check if scan already running %
set(0,'showhiddenhandles', 'on');
hAllFigs = get(0, 'children');
hFig = findobj(hAllFigs, 'flat', 'tag', 'scanBaseWindow');
set(0,'showhiddenhandles', 'off');
if ~isempty(hFig);
    figure(hFig);   return
end

% App data struct: initialise these fields %
scanData.selData = {};
scanData.allData = {};
scanData.selTrial = [];
scanData.selCell = [];
scanData.selMap = [];
scanData.infoDataShow = [];
scanData.infoDataSort = [];
scanData.infoTrialShow = [];
scanData.infoTrialSort = [];
scanData.infoCellShow = [];
scanData.infoCellSort = [];
% App data struct: load these fields in already %
app_path = getenv('appdata');
try
    L = load([app_path, '\matlab_scan', '\settings.mat']);
    scanData.settings = L.settings;
catch
    % If failed, make a default structure %
    warning('SCAn:base_window', 'Couldn''t load saved settings');
    scanData.settings.lastDirLoad = pwd;
    scanData.settings.lastDirOpen = pwd;
    scanData.settings.lastDirSave = pwd;
    scanData.settings.lastDirEps = pwd;
    scanData.settings.scanWindowFormat = 'normal';     
    scanData.settings.fixedPlotHeight = [];
    scanData.settings.maxRowNumber = [];
    scanData.settings.defaultCutTag = '';
end
% The following are fields that have been added in later versions: need to check for existence, then assign if not there %
if ~isfield(scanData.settings,'scanWindowFormat');    scanData.settings.scanWindowFormat = 'normal';    end
if ~isfield(scanData.settings,'fixedPlotHeight');     scanData.settings.fixedPlotHeight = [];           end
if ~isfield(scanData.settings,'maxRowNumber');        scanData.settings.maxRowNumber = [];              end
if ~isfield(scanData.settings,'defaultCutTag');       scanData.settings.defaultCutTag = '';             end

% Look for format argument at input line %
if ~isempty(varargin)
    if strcmp(varargin{1},'normal')
        scanData.settings.scanWindowFormat = 'normal';
    elseif strcmp(varargin{1},'tall')
        scanData.settings.scanWindowFormat = 'tall';
    end
end     
% App data struct: Get screen size (Get every time) %
reset_str = get(0, 'units');
set(0, 'units', 'characters');
scanData.screenChar = get(0, 'screensize');
set(0, 'units', 'pixels');
scanData.screenPix = get(0, 'screensize');
set(0, 'units', reset_str);

% Open figure %
if strcmp(scanData.settings.scanWindowFormat,'tall')
    hFig = openfig('scan_base_tall_format.fig');
else
    hFig = openfig('scan_base.fig');
end
H = guihandles(hFig);

% Position Figure %
figPos = get(hFig,'position');
set(hFig, 'position', [(scanData.screenChar(3:4)-figPos(3:4))./2 figPos(3:4)]);

% App data struct: save figure handle %
scanData.figHandle = hFig;
scanData.same_cells = 0;    % This is for backwards compatibility.
% App data struct: save %
guidata(hFig, scanData);

% Put extra menus onto base window %
scan_base('baseMenus',hFig);

% Populate list boxes with existing data %
scan_base('listRefreshData', H.listData);