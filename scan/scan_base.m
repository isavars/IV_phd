function [varargout] = scan_base(funStr, varargin)
%       scan_base('callbackSubfunction', CallingObjectHandle);
%       scan_base('callbackSubfunction');
%
% Can also output a return argument:
%
%       rtn = scan_base('callbackSubfunction');
%
% Contains the callbacks of most of the components of the base window:
% those that deal with GUI, file IO and data struct manipulation.
%
% Callbacks relating to actual data processing have their own functions.

% Contents (functions):
%   baseMenus
%   openData
%   saveData
%   closeData
%   closeApp
%   listRefreshData
%   listRefreshMap
%   listRefreshTrial
%   listRefreshCell
%   displayInfo
%   drawPreviewMap
%   ratesDataHandler
%   ratesMapNamer
%   structRemoveMaps
%   structRemoveTrialCell
%   structAddHistology
%   saveEps
%   savePdf
%   optionsDialog
%   userVar
%   userVarRedraw

%%% Call relevant subfunction %%%
if ~isempty(varargin)
    hObj = varargin{1};
else
    SD = gss;   hObj = SD.figHandle;  % If no hObj passed, just get base figure handle.
end
hFun = str2func(funStr);
if nargout==0
    hFun(hObj, guihandles(hObj), guidata(hObj));
elseif nargout==1
    varargout{1} = hFun(hObj, guihandles(hObj), guidata(hObj));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = baseMenus(hObj, H, SD)
%%% Put menus on the base figure. Coding them in here is easier then using GUIDE. %%% 
%%% 'Dataset' Menu %%%
uimenu('label','Load Histology Image','parent',H.menuData,'callback','scan_base(''structAddHistology'');');
%%% 'Figures' Menu %%%
hMenuFig = uimenu('label', 'Figures', 'parent', hObj);
    uimenu('label', 'Draw all cell properties', 'parent', hMenuFig, 'callback', 'gra_allplots(gss);');
    hTempAC=uimenu('label', 'Autocorrelogram', 'parent', hMenuFig);
        uimenu('label', '20ms', 'parent', hTempAC, 'callback', 'gra_xcfig(gss,''xcBin'',0.001,''xcWin'',0.02);');
        uimenu('label', '500ms', 'parent', hTempAC, 'callback', 'gra_xcfig(gss,''xcBin'',0.01,''xcWin'',0.5);');
    hTempCC=uimenu('label', 'Cross-correlogram', 'parent', hMenuFig);
        uimenu('label', '20ms', 'parent', hTempCC, 'callback', 'gra_xcfig(gss,''corrMode'',''cross'',''xcBin'',0.001,''xcWin'',0.02);');
        uimenu('label', '50ms', 'parent', hTempCC, 'callback', 'gra_xcfig(gss,''corrMode'',''cross'',''xcBin'',0.001,''xcWin'',0.05);');
        uimenu('label', '100ms', 'parent', hTempCC, 'callback', 'gra_xcfig(gss,''corrMode'',''cross'',''xcBin'',0.001,''xcWin'',0.1);');
        uimenu('label', '500ms', 'parent', hTempCC, 'callback', 'gra_xcfig(gss,''corrMode'',''cross'',''xcBin'',0.001,''xcWin'',0.5);');
        uimenu('label', 'Cross-corr CL 2013', 'parent', hTempCC, 'callback', ...
            ['gra_acfig(gss,''corrMode'',''cross'',''spikingProbYAxis'',1,''acWin'',50,''acBin'',1,''sigThr'',''3SD'');' ...
             'gra_acfig(gss,''corrMode'',''cross'',''spikingProbYAxis'',1,''acWin'',50,''acBin'',1,''sigThr'',''3SD'',''xCorrSubtraction'',''multiShift'');']  );
    hGridProps=uimenu('label', 'Grid Props + AutoCorrs', 'parent', hMenuFig);
        uimenu('label', 'Sargolini 2006', 'parent', hGridProps, 'callback', 'gra_gridprops(gss,''Sargolini 2006'');');
        uimenu('label', 'Barry 2007', 'parent', hGridProps, 'callback', 'gra_gridprops(gss,''Barry 2007'');');
        uimenu('label', 'Wills 2010', 'parent', hGridProps, 'callback', 'gra_gridprops(gss,''Wills 2010'');');
        uimenu('label', 'Lever 2012', 'parent', hGridProps, 'callback', 'gra_gridprops(gss,''Lever 2012'');');
        uimenu('label', 'Watershed', 'parent', hGridProps, 'callback', 'gra_gridprops(gss,''Watershed'');');
    uimenu('label', 'Time series maps', 'parent', hMenuFig, 'callback', 'gra_mapfig_series;');    
    uimenu('label', '2-D FFT', 'parent', hMenuFig, 'callback', 'gra_fft2(gss,[]);');
    uimenu('label', 'Directional Distributive Hypothesis', 'parent', hMenuFig, 'callback', 'gra_distributivehypothesis(gss);');
    uimenu('label','Spatial Information','parent',hMenuFig,'callback','gra_spatialinfo(gss)');
    uimenu('label', 'EEG', 'parent', hMenuFig, 'callback', 'gra_eeg(gss);');
    hMenuFigWf = uimenu('label', 'Plot Waveforms', 'parent', hMenuFig);
        uimenu('label', 'Max Channel', 'parent', hMenuFigWf, 'callback', 'gra_wf(''maxchannel'');');
        uimenu('label', 'All Channels', 'parent', hMenuFigWf, 'callback', 'gra_wf(''allchannels'');');
    uimenu('label','Time Window AC','parent',hMenuFig,'callback','gra_timewindowdisphist(gss);');
    uimenu('label', 'Dir Props', 'parent', hMenuFig, 'callback', 'gra_dirprops();');
%%% Analyse Data %%%
hMenuAnal =  uimenu('label', 'Analyse', 'parent', hObj);
    uimenu('label','Cluster Isolation','parent',hMenuAnal,'callback','analysis_clusterisolation(gss);');
    uimenu('label','Output to Excel','parent',hMenuAnal,'callback','analysis_output2excel(gss);');
    % Developmental project %
    hMenuExpDev = uimenu('label','Development','parent',hMenuAnal);
        uimenu('label','Plot Speed-matched FFT','parent',hMenuExpDev,'callback','dev_fftplot(gss);');
        uimenu('label','Theta Mod','parent',hMenuExpDev,'callback','spk_thetamod(gss);');
        
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = openData(hObj, H, SD)
%%% Open data variables into workspace %%%
dir_reset = pwd;
try   cd(SD.settings.lastDirOpen);
catch;   SD.settings.lastDirOpen = pwd;
end
[filenames, pathname] = uigetfile({'*.mat'}, 'Open Data .. ', 'multiselect', 'on');
if ~ischar(pathname)
    cd(dir_reset);   return
end
filenames = cellstr(filenames);
for ii=1:length(filenames)   
	S = load([pathname, filenames{ii}]);
    names = fieldnames(S);    
    for jj=1:length(names)        
        assignin('base', names{jj}, S.(names{jj}));        
    end
end
SD.settings.lastDirOpen = pathname;
cd(dir_reset);
guidata(hObj,SD);
scan_base('listRefreshData', H.listData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = saveData(hObj, H, SD)
% Save variable, prompting for filename %
dir_reset = pwd;
try   cd(SD.settings.lastDirSave);
catch;   SD.settings.lastDirSave = pwd;
end
if length(SD.selData) == 1
    fname = SD.selData{1};
else
    fname = 'scan';
end
[filename, pathname] = uiputfile([fname, '.mat'], 'Save Dataset as ..');
cd(dir_reset);
if filename==0;  return;  end;
SD.settings.lastDirSave = pathname;
% Save. This is in format SAVE('filename','var1','var2'), as windows desktop path makes trouble otherwise.
% Watch out for the string arguments: <''''> = <'> in a string, <'''> = <'> at the end of a string.
dnamer = [', ''', SD.selData{1}, ''''];
if length(SD.selData) > 1
    for ii=2:length(SD.selData)
        dnamer = [dnamer, ', ''', SD.selData{ii}, ''''];
    end
end
evalin('base', ['save(''', [pathname, filename], '''', dnamer, ');']);
guidata(hObj,SD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = closeData(hObj, H, SD)
%%% Clear data set(s) from workspace (NB. This is now called 'Delete' on menu) %%%
if ~isempty(SD.selData)
    dnamer = SD.selData{1};
    for ii=2:length(SD.selData)
        dnamer = [dnamer, '   ', SD.selData{ii}];
    end
    check = questdlg(['Are you sure you want to delete dataset(s):  ', dnamer], 'SCAn: Delete datasets', 'OK', 'Cancel', 'OK');
    if strcmp(check, 'OK')
        evalin('base', ['clear ',  dnamer, ';']);
    end
    scan_base('listRefreshData', H.listData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = closeApp(hObj, H, SD)
%%% Close SCAn app - save settings from SCAn structure %%%
app_path = getenv('appdata');
status = 1;
if isempty(dir([app_path '\matlab_scan']))
    status = mkdir(app_path, 'matlab_scan');  
end
if status == 1
    settings = SD.settings;
    save([app_path, '\matlab_scan\settings.mat'], 'settings');
else
    disp('SCAn: Didn''t save settings. Couldn''t create directory.');
end
delete(SD.figHandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = listRefreshData(hObj, H, SD)
%%% Refresh data list box on base window %%%
% This function also starts a cascade to refresh all boxes:
% mapListRefresh -> trialListRefresh -> cellListRefresh -> drawPreviewMap
all_list = evalin('base', 'who');   dataList = {};
for ii=1:length(all_list)
    if evalin('base', ['isfield(', all_list{ii}, ',''is_scan'')'])
        dataList{end+1} = all_list{ii};
    end
end
hObj = H.listData;  % This function is used by external calling functions, so hObj input might be base figure handle.
SD.allData = dataList;
if isempty(dataList) 
    % No SCAN structs in workspace - set list blank %
    set(hObj,'string',{''},'value',1);   
    set(H.infoData,'string',{''},'value',1);
    SD.allData = [];
else
    % If there is data .. %
    % Compile display list with info var, if necessary %
    if isempty(SD.infoDataShow)
        dataListDisp = dataList;
    else
        for ii=1:length(dataList)
            if evalin(   'base',  ['isfield(' dataList{ii} '.user, ''' SD.infoDataShow ''')']   )
                dataListDisp{ii} = [dataList{ii}, '   ::   ', evalin('base', [dataList{ii} '.user.' SD.infoDataShow])];
            else
                dataListDisp{ii} = [dataList{ii}, '   ::'];
            end
        end
    end
    % Sort %
    % NOTE, 12/11/08. This is a quick and dirty sort implementation. 
    % 1) Doesn't worry about maintaining selection focus during sort.
    % 2) For data (not for trials or cells) don't need to worry about
    % matching display index to 'real' data index (as data selection based
    % on string matching).
    if ~isempty(SD.infoDataSort)
        for ii=1:length(dataList)
            if evalin( 'base',  ['isfield(' dataList{ii} '.user, ''' SD.infoDataShow ''')']  )
                sortList{ii} = evalin('base', [dataList{ii} '.user.' SD.infoDataShow]);
            else
                sortList{ii} = '';
            end
        end
        [dummy sortInd] = sort(sortList);
        dataListDisp = dataListDisp(sortInd);
        dataList = dataList(sortInd);
    end
    % overflow check %
    if length(dataList) < max(get(hObj, 'value')) 
        set(hObj, 'value', length(dataList));     
    end                                            
    set(hObj, 'string', dataListDisp);             
    SD.selData = dataList(get(hObj, 'value'));
    % .. and Data Info %
    if length(SD.selData)==1
        str{1} = [SD.selData{1} ':'];
        data = evalin('base', SD.selData{1});   f = fieldnames(data.user); 
        for ii=1:length(f);  
            if ~isempty(data.user.(f{ii})) && ischar(data.user.(f{ii}));   str{end+1} = [f{ii} ' = ' data.user.(f{ii})];   end
        end
    else   str = {'Multiple Data'};
    end
    if length(str) < max(get(H.infoData, 'value'))  % 
        set(H.infoData, 'value', length(str));      % overflow check
    end                                             %  
    set(H.infoData, 'string', str);
end
% Refresh Other List Boxes %
guidata(hObj, SD); 
scan_base('listRefreshMap', H.listMap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = listRefreshMap(hObj, H, SD)
%%% Refresh map list box on base window %%%
% selMap is not an index integer but, rather, a .params rate map struct
% Need to go through all selected datasets, a create a (unique) list of 
% maps exist - mapList.This is a metadata-only .rate_maps struct.
if isempty(SD.selData)
    set(hObj,'string',{''},'value',1);   SD.selMap = [];
elseif length(SD.selData) > 1
    set(hObj,'string',{'Multiple datasets selected'},'value',1);   SD.selMap = [];
else
    % Make unique list of maps in all selected datasets %
    mapList = [];
    data=evalin('base', SD.selData{1});
    for ii=1:length(data.rate_maps)
    	mapList(ii).name = data.rate_maps(ii).name;
    end
    % Set list box and map selection %
    if isempty(mapList)
        set(hObj, 'string', '', 'value', 1);   SD.selMap = [];
    else
        [mapListDisp{1:length(mapList)}] = deal(mapList.name);
        if max(get(hObj, 'value')) > length(mapListDisp) % 
            set(hObj, 'value', length(mapListDisp));     % Overflow check
        end                                              %
        set(hObj, 'string', mapListDisp);
        SD.selMap = get(hObj, 'value');   
    end
end
% Refresh Cell Box %
guidata(hObj, SD);
scan_base('listRefreshTrial', H.listTrial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = listRefreshTrial(hObj, H, SD)
%%% Refresh map list box on base window %%%
% Also calls cellListRefresh, to make infoCell box relate to selected trial.
if isempty(SD.selData)
    % No data in workspace - all blank %
    set(hObj, 'string', {''}, 'value', 1);
    set(H.infoTrial, 'string', {''}, 'value', 1);
    SD.selTrial = [];
elseif length(SD.selData)>1
    % If more than one data selected - don't list anything %
    set(hObj,'string','Multiple Data Selected', 'value', 1);   
    set(H.infoTrial, 'string', 'Multiple Data Selected', 'value', 1);   
    SD.selTrial=[];
elseif length(SD.selData)==1
    % One data struct selected: list trials %
    data = evalin('base', SD.selData{1});
    [trialList{1:length(data.trials)}] = deal(data.trials.trialname);
    % Compile display list with info var, if necessary %
    if isempty(SD.infoTrialShow)
        trialListDisp = trialList;
    else
        for ii=1:length(trialList)
            switch SD.infoTrialShow
                case 'trial_length'
                    trialListDisp{ii} = [trialList{ii} '   ::   ' num2str(data.trials(ii).dur)];
                case 'original_ppm'
                    trialListDisp{ii} = [trialList{ii} '   ::   ' num2str(data.trials(ii).original_ppm)];
                otherwise
                    if isfield(data.trials(ii).user,SD.infoTrialShow)
                        trialListDisp{ii} = [trialList{ii} '   ::   ' data.trials(ii).user.(SD.infoTrialShow)];
                    else
                        trialListDisp{ii} = [trialList{ii} '   ::'];
                    end
            end
        end             
    end
    trialListDisp = ['<Select All Trials>', trialListDisp];
    if max(get(hObj, 'value')) > length(trialListDisp)  % 
        set(hObj, 'value', length(trialListDisp));      % check for value overflow 
    end                                                 %
    set(H.listTrial, 'string', trialListDisp);
    if get(hObj, 'value') == 1;
        SD.selTrial = (1:length(data.trials)); % For <Select All Trials> 
    else
        SD.selTrial = get(hObj, 'value') - 1;  % For 1 or more trials selected
    end
    % Trial Info %
    if length(SD.selTrial) == 1
        str{1} = [data.trials(SD.selTrial).trialname ':'];
        if isstruct(data.trials(SD.selTrial).eeg)
            if isempty(data.trials(SD.selTrial).eeg(1).eeg)
                str{end+1} = 'No EEG Found!';
            else
                for ii=1:length(data.trials(SD.selTrial).eeg)
                    eeg = data.trials(SD.selTrial).eeg(ii);
                    filtStr = {'', 'Direct+Notch', '', '', 'LP+Notch',''};
                    str{end+1} = ['EEG' num2str(ii) ': Ch ' num2str(eeg.channel) ', ' num2str(round(eeg.scalemax)) 'uV, ' filtStr{eeg.filter+1}];
                end
            end
        end
        str{end+1} = ['trial_length = ' num2str(data.trials(SD.selTrial).dur)];
        str{end+1} = ['ppm = ' num2str(data.trials(SD.selTrial).ppm)];
        str{end+1} = ['original_ppm = ' num2str(data.trials(SD.selTrial).original_ppm)];
        f = fieldnames(data.trials(SD.selTrial).user);
        for ii=1:length(f);  
            if ~isempty(data.trials(SD.selTrial).user.(f{ii})) && ischar(data.trials(SD.selTrial).user.(f{ii}));
                str{end+1} = [f{ii} ' = ' data.trials(SD.selTrial).user.(f{ii})];
            end
        end    
    else
        str = {'Multiple Trials Selected'};
    end
    if length(str) < max(get(H.infoTrial, 'value'))  % 
        set(H.infoTrial, 'value', length(str));      % overflow check
    end                                             %  
    set(H.infoTrial, 'string', str);
end
% Refresh Cell Box %
guidata(hObj, SD);   
scan_base('listRefreshCell', H.listCell);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = listRefreshCell(hObj, H, SD)
%%% Refresh map list box on base window %%%
if isempty(SD.selData)
    % No data in workspace - all blank %
    set(hObj,'string',{''}, 'value', 1);   
    set(H.infoCell, 'string', {''}, 'value', 1);   
    SD.selTrial=[];   return
end
if length(SD.selData)>1
    % If more than one data selected, don't list anything %
    set(hObj,'string','Multiple Data Selected', 'value', 1);   
    set(H.infoCell, 'string', 'Multiple Data Selected', 'value', 1);   
    SD.selTrial=[];   return
end
data = evalin('base', SD.selData{1});
if isempty(data.trials(1).cells)
    set(hObj, 'string', {'No Cells'}, 'value', 1);
    set(H.infoCell, 'string', {'No Cells'}, 'value', 1);
    SD.selCell = [];
elseif ~isempty(data.rate_maps) && strcmp(data.rate_maps(SD.selMap).params.mode, 'pos');
    % Pos maps - special case %
    set(hObj, 'string', {'Position Maps'}, 'value', 1);
    set(H.infoCell, 'string', {'Position Maps'}, 'value', 1);
    SD.selCell = 1;
else
    % Normal cell list %
    for ii=1:length(data.trials(1).cells)
        cellList{ii} = [num2str(ii) ': Cell ', num2str(data.trials(1).cells(ii).cellnum), ' Tet ', num2str(data.trials(1).cells(ii).tet)];
    end
    % Compile display list with info var, if necessary %
    if isempty(SD.infoCellShow)
        cellListDisp = cellList;
    else
        for ii=1:length(cellList)
            if strcmp(SD.infoCellShow, 'N_Spike')
                varStr = num2str(length(data.trials(SD.selTrial).cells(ii).st));
            elseif ~ischar(data.trials(1).cells(ii).user.(SD.infoCellShow))
                varStr = num2str(data.trials(1).cells(ii).user.(SD.infoCellShow));
            else
                varStr = data.trials(1).cells(ii).user.(SD.infoCellShow);
            end
            if length(SD.selTrial)==1
                cellListDisp{ii} = [cellList{ii} '   ::   ' varStr];
            else
                cellListDisp{ii} = [cellList{ii} '   ::   '];
            end
        end             
    end
    cellListDisp = ['<Select All Cells>' cellListDisp]; 
    if max(get(hObj, 'value')) > length(cellListDisp) % 
        set(hObj, 'value', length(cellListDisp));     % overflow check
    end                                               %
    set(hObj, 'string', cellListDisp);           
    if get(hObj, 'value') == 1;
        SD.selCell = (1:length(data.trials(1).cells)); % For <Select All Cells>
    else
        SD.selCell = get(H.listCell, 'value') - 1;     % For 1 or more cells selected
    end
    % Cell Info %
    if length(SD.selCell)==1 && length(SD.selTrial)==1
        str{1} = ['Cell ', num2str(data.trials(1).cells(SD.selCell).cellnum), ' Tet ', num2str(data.trials(1).cells(SD.selCell).tet) ':'];
        str{2} = ['N_Spike = ' num2str(length(data.trials(SD.selTrial).cells(SD.selCell).st))];
        str{3} = ['Cut = ' data.trials(SD.selTrial).cells(SD.selCell).cut];
        % Caution: some old structs have cells.user=[], not struct([])
        try  f = fieldnames(data.trials(SD.selTrial).cells(SD.selCell).user);  catch;  f = {};  end
        for ii=1:length(f);
            if ~isempty(data.trials(SD.selTrial).cells(SD.selCell).user.(f{ii}))
                if ischar(data.trials(SD.selTrial).cells(SD.selCell).user.(f{ii}))
                    str{end+1} = [f{ii} ' = ' data.trials(SD.selTrial).cells(SD.selCell).user.(f{ii})];
                elseif isnumeric(data.trials(SD.selTrial).cells(SD.selCell).user.(f{ii}))
                    str{end+1} = [f{ii} ' = ' num2str(data.trials(SD.selTrial).cells(SD.selCell).user.(f{ii}))];
                else
                    continue
                end
            end
        end
    else
        str = {'Multiple Trials/Cells'};
    end
    if length(str) < max(get(H.infoCell, 'value'))  % 
        set(H.infoCell, 'value', length(str));      % overflow check
    end  
    set(H.infoCell, 'string', str);
end
guidata(hObj,SD);
% Draw into preview axes %
scan_base('drawPreviewMap', H.previewAxes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = displayInfo(hObj, H, SD)
%%% Add 'info' metadata tags to the data, trial and cell lists.  %%%
% Ths function is only a) sets the show/sort fields in ScanData, b) alters
% the list title box. Then calls listRefreshXxx, which actually modifies 
% list box contents.
tag = get(hObj,'tag');
hInfo = H.(['info' tag(5:end)]);
if get(hInfo,'value')==1
    infoVar = [];            % Top in list is object label. Doing 'Show' 
    titleStr = tag(5:end);   % on this will remove info from list box.
else
    infoList = get(hInfo,'string');
    infoVar = strtok(infoList{get(hInfo,'value')});
    titleStr = [tag(5:end) '   ::   ' infoVar];
end
SD.(['info' tag(5:end) 'S' tag(2:4)]) = infoVar;  % Set a field in ScanData, with format (eg) SD.infoDataSort 
set(H.(['listTitle' tag(5:end)]), 'string', titleStr);      % Set list title
% Refresh Lists %
guidata(hObj,SD);   
scan_base(['listRefresh' tag(5:end)], H.(['list' tag(5:end)]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = drawPreviewMap(hObj, H, SD)
%%% Draw a map into the preview box (or text, if no single map selected) %%
delete(get(hObj, 'children'));
if get(H.previewMap,'value')
    if isempty(SD.selCell) ||isempty(SD.selMap) || length(SD.selTrial)>1 || length(SD.selCell)>1 || length(SD.selCell)>1
        text('string', {'Select a trial and cell',' to show a map'}, 'position', [0.5 0.5], 'units', 'normalized', ...
            'horizontalalignment', 'center', 'parent', hObj);
    else
        data = evalin('base', SD.selData{1});
        mInd = SD.selMap;
        gra_plotmap(data.rate_maps(mInd).maps{SD.selTrial, SD.selCell}, 'handle', hObj);
    end
elseif get(H.previewHistology,'value')
    if length(SD.selData) == 1
        data = evalin('base', SD.selData{1});
        if isfield(data,'histology') && ~isempty(data.histology)
            image(data.histology,'parent',hObj);
            axis(hObj,'off');
            set(hObj,'tag','previewAxes');  % As with rate maps, image command wipes out axes tag.
        else
            text('string', {'Use Dataset > Load Histology',' to load a histology image'}, 'position', [0.5 0.5], 'units', 'normalized', ...
            'horizontalalignment', 'center', 'parent', hObj);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = ratesDataHandler(hObj, H, SD)
%%% Interface between GUI and RATES_MAIN. Calls RATES_UIPARAMS, passes
% output and dataset struct to RATES_MAIN, assigns resulting maps cell array, 
% params struct and name string to data.rate_maps(end+1). Useful for:
% - when making the same rate maps in multiple datasets
% - when making map series (each params(n) passed individually in a loop)
% - when transforming maps (first make maps, then call TRANS_MAIN) 
%%% TODO - trnasforming maps - copy in form old code.
params=rates_params;
if isempty(params);   return;   end     % Cancel button pressed

%%%%%%%%%% INSERT HACKS FOR ODD MAPS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params.filt_index = [1:300*50 600*50:900*50];

% if params.smooth==5
%     params.smooth=[0.0025 0.0125 0.0200 0.0125 0.0025;...
%                    0.0125 0.0625 0.1000 0.0625 0.0125;...
%                    0.0200 0.1000 0.1600 0.1000 0.0200;...
%                    0.0125 0.0625 0.1000 0.0625 0.0125;...
%                    0.0025 0.0125 0.0200 0.0125 0.0025;];
% end
% disp(size(params.smooth));
                   
for ii=1:length(SD.selData)
    data = evalin('base', SD.selData{ii});
    nextMapInd=length(data.rate_maps)+1;
    if strcmp(params.alg,'pxd')
        [data.rate_maps(nextMapInd).maps data.rate_maps(nextMapInd+1).maps] = rates_main(data, params);
        dParams=params; dParams.space='dir';                            % In PxD params, space='place'. For consistency,
        data.rate_maps(nextMapInd+1).name = ratesMapNamer(dParams);     % set to 'dir' for the dir maps (2nd in
        data.rate_maps(nextMapInd+1).params = dParams;                  % the struct).
    elseif ~isempty(params.adaptive_smooth) && strcmp(params.mode,'pos');
        [dum,data.rate_maps(nextMapInd).maps] = rates_main(data,params);
    else
        data.rate_maps(nextMapInd).maps = rates_main(data,params);
    end
    % Transform maps, if required %
    if params.trans_active
        % Check if shape codes have been assigned to all trials %
        for jj=1:length(data.trials);   
            if isempty(data.trials(jj).user.shape);  error('Can''t transform maps - you need to set the ''shape'' variable (8=square, 1=circle, etc ) for all trials');   end;  
        end
        % Transform maps %
        maps=data.rate_maps(nextMapInd).maps;
        transMaps=cell(size(maps));
        targetMap=maps{params.trans_target,1};
        targetShape=str2double(data.trials(  params.trans_target  ).user.shape);
        for jj=1:length(data.trials)
            transMaps(jj,:) = rates_transform(maps(jj,:),str2double(data.trials(jj).user.shape),targetMap,targetShape,'reg_size',params.trans_reg_size);
        end
        data.rate_maps(nextMapInd).maps = transMaps;
    end
    % Assign map metadata in data struct, and then in data base workspace %
    data.rate_maps(nextMapInd).name = ratesMapNamer(params);
    data.rate_maps(nextMapInd).params = params;
    assignin('base', SD.selData{ii}, data);
end
% If called with just one data struct selected, draw the maps %
if length(SD.selData)==1
    maps = data.rate_maps(nextMapInd).maps;
    gra_mapfig(data,nextMapInd,1:size(maps,1),1:size(maps,2));
end
scan_base('listRefreshData', H.listData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [str] = ratesMapNamer(params)
%%% Generate a name string for a set of rate maps %%%
% Check for map series, generate name if needed %
str = '';
filt_check = 'no_filts';
if length(params) > 1
    filts = {'time', 'speed', 'dir'};
    for ii=1:length(filts)
        f1 = params(1).(['filt_' filts{ii}]);
        f2 = params(2).(['filt_' filts{ii}]);
        if ~isempty(f1)
            if f1~=f2
                f3 = params(end).(['filt_' filts{ii}]);
                str = [filts{ii}, ' series: [', num2str(f1(1)), ':', num2str(f1(2)-f1(1)), ':', num2str(f3(2)), '], '];
                filt_check = filts{ii};
                break
            end
        end
    end
end
% Now all the other parameters: only give name if value not default %
str = [str, params(1).space, '-', params(1).mode];
if ~strcmp(params(1).alg, 'pd')
    str = [str, '-' params(1).alg];
end
if strcmp(params.space,'place')     % Always put the bin size
    str = [str, ', bin=', num2str(params(1).bin)];
elseif strcmp(params.space,'dir')
    str = [str, ', bin=', num2str(params(1).bin_dir)];
end
if strcmp(params.space,'place') && (max(size(params.smooth))>1 || params(1).smooth~=5  || ~isempty(params.adaptive_smooth))
    if isempty(params.adaptive_smooth)
        str = [str, ', smooth=', num2str(params(1).smooth)];
    else
        str = [str, ', smooth=adapt:' num2str(params(1).adaptive_smooth)];
    end
elseif strcmp(params.space,'dir') && params(1).smooth_dir~=5
    str = [str, ', smooth=', num2str(params(1).smooth_dir)];
end
% Don't put filter params when they are part of a series %
if ~strcmp(filt_check, 'time')
	if ~isempty(params(1).filt_time)
        str = [str, ', time=', num2str(params(1).filt_time(1,1)), '-', num2str(params(1).filt_time(1,2)), 's'];
	end
end
if ~strcmp(filt_check, 'speed')
	if ~isempty(params(1).filt_speed)
        str = [str, ', speed=', num2str(params(1).filt_speed(1)), '-', num2str(params(1).filt_speed(2)), 'cm\s'];
	end
end
if ~strcmp(filt_check, 'dir')
	if ~isempty(params(1).filt_dir)
        str = [str, ', dir=', num2str(params(1).filt_dir(1)), '-', num2str(params(1).filt_dir(2)), 'deg'];
	end
end
if ~strcmp(filt_check, 'x')
	if ~isempty(params(1).filt_x)
        str = [str, ', x=', num2str(params(1).filt_x(1)), '-', num2str(params(1).filt_x(2))];
	end
end
if ~strcmp(filt_check, 'y')
	if ~isempty(params(1).filt_y)
        str = [str, ', y=', num2str(params(1).filt_y(1)), '-', num2str(params(1).filt_y(2))];
	end
end
% Crop %
if ~isempty(params.crop) && strcmp(params.space,'place')
    str = [str, ', crop=', num2str(params.crop(1)), ',' num2str(params.crop(2))];
end
% Transform %
if params.trans_active
    str=[str ', trans_to_tr' num2str(params.trans_target), ', sameSize=' ];
    if params.trans_reg_size;   str=[str 'yes'];    else   str=[str 'no'];    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = structRemoveMaps(hObj, H, SD)
%%% Delete rate maps from data struct %%%
userCheck = questdlg('Are you sure want to delete the selected maps?', 'Delete Maps', 'Yes, delete', 'Cancel', 'Yes, delete');
if strcmp(userCheck, 'Yes, delete')
    for ii=1:length(SD.selData)
        data = evalin('base', SD.selData{ii});
        if ~isempty(data.rate_maps)
            rmInd = SD.selMap;
            if ~isempty(rmInd)
                data.rate_maps = data.rate_maps( setdiff(1:length(data.rate_maps),rmInd) );
            end
        end
        assignin('base', SD.selData{ii}, data);
    end
end
scan_base('listRefreshData', H.listData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = structRemoveTrialCell(hObj, H, SD)
%%% Remove trials or cells from selected dataset %%%
if length(SD.selData)>1
    errordlg({'To remove Trials/Cells, please select only one dataset'},'SCAn Error');
end
data = evalin('base', SD.selData{1});
hList = get(SD.figHandle, 'currentobject');   tagList = get(hList,'tag');
switch tagList(5:end)
    case 'Trial'
        namer = '';
        for ii=1:length(SD.selTrial)
            namer = [namer data.trials(SD.selTrial(ii)).trialname '  '];
        end
        userCheck = questdlg(['Definitely remove trials ' namer '?'], 'Remove Trials', 'Yes, remove', 'Cancel', 'Yes, remove');
        if strcmp(userCheck, 'Yes, remove')
            rmInd = setdiff(1:length(data.trials),SD.selTrial);
            data.trials = data.trials(rmInd);
            for ii=1:length(data.rate_maps)
                data.rate_maps(ii).maps = data.rate_maps(ii).maps(rmInd,:);
            end
        end
    case 'Cell'
        namer = '';
        for ii=1:length(SD.selCell)
            c = data.trials(1).cells(SD.selCell(ii));
            namer = [namer 'Cell' num2str(c.cellnum) ' Tet' num2str(c.tet) '   '];
        end
        userCheck = questdlg(['Definitely remove cells  ' namer '?'], 'Remove Cells', 'Yes, remove', 'Cancel', 'Yes, remove');
        if strcmp(userCheck, 'Yes, remove')
            rmInd = setdiff(1:length(data.trials(1).cells),SD.selCell);
            for ii=1:length(data.trials)
                data.trials(ii).cells = data.trials(ii).cells(rmInd);
            end
            for ii=1:length(data.rate_maps)
                if ~strcmp(data.rate_maps(ii).params.mode,'pos')
                    data.rate_maps(ii).maps = data.rate_maps(ii).maps(:,rmInd);
                end
            end
        end
end
assignin('base', SD.selData{1}, data);
scan_base('listRefreshData',H.listData);      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = structAddHistology(hObj,H,SD)
% Load an image and store in data.histology field. Will downsample image to
% 275 x 275 pixels (size of preview window). Also set to int8?
dir_reset = pwd;
try
    cd(SD.settings.lastDirHisto);
catch
    SD.settings.lastDirHisto = pwd;
end
formatList={'*.jpg;*.jpeg', 'jpg'; '*.tif;*.tiff', 'tiff'};
[filename, pathname formatIndex] = uigetfile(formatList, 'Select Histology Image .. ', 'multiselect', 'off');
if ~ischar(pathname)
    cd(dir_reset);   return
end
Im = imread([pathname, filename],formatList{formatIndex,2});
for ii=1:length(SD.selData)
    data=evalin('base',SD.selData{ii});
    data.histology = Im;
    assignin('base',SD.selData{ii},data);
end
SD.settings.lastDirHisto = pathname;
guidata(SD.figHandle, SD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = saveEps(hObj,H,SD)
% Save a figure as an eps, prompting for filename.
% Note that hObj is the handle of the figure to print - need to get SD directly.
SD = gss;
dir_reset = pwd;
figure(hObj);
try
    cd(SD.settings.lastDirEps);
catch
    SD.settings.lastDirEps = pwd;
end
% Figure title is file name - need to deblank %
def_title = get(gcf, 'name');
for ii=1:length(def_title)
    if strcmp(def_title(ii), ' ')
        def_title(ii) = '_';
    end
end
cd(SD.settings.lastDirEps);
[filename, pathname] = uiputfile( [def_title, '.eps'], 'Save figure as eps ..');
cd(dir_reset);
if filename==0
    return
end
% This is in format PRINT('filename','var1','var2') because the windows desktop path makes trouble otherwise.
set(hObj,'renderer','painters');   % To make sure all figures outputtd as vector graphics
evalin('base', [   'print(gcf, ''-depsc'', ''', [pathname, filename], ''');'    ]     );
if ~isempty(SD)
    SD.settings.lastDirEps = pathname;
    guidata(SD.figHandle, SD);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = savePdf(hObj,H,SD)
% Save a figure as a pdf, prompting for filename.
% Note that hObj is the handle of the figure to print - need to get SD directly.
SD = gss;
dir_reset = pwd;
figure(hObj);
try
    cd(SD.settings.lastDirEps);
catch
    SD.settings.lastDirEps = pwd;
end
% Figure title is file name - need to deblank %
def_title = get(gcf, 'name');
for ii=1:length(def_title)
    if strcmp(def_title(ii), ' ')
        def_title(ii) = '_';
    end
end
cd(SD.settings.lastDirEps);
[filename, pathname] = uiputfile( [def_title, '.pdf'], 'Save figure as pdf ..');
cd(dir_reset);
if filename==0
    return
end

set(hObj,'renderer','painters');   % To make sure all figures outputtd as vector graphics
% Make sure the figure fits on the page.
% sz = get(hObj,'paperposition');
% if sz(3) > hObj.PaperSize(1)
%     sz([3 4]) = sz([3 4]) .* (sz(3) / hObj.PaperSize(1));
% end
% if sz(4) > hObj.PaperSize(2)
%     sz([3 4]) = sz([3 4]) .* (sz(4) / hObj.PaperSize(4));
% end
% sz([1 2]) = 0;
% set(hObj,'paperposition',sz);
% This is in format PRINT('filename','var1','var2') because the windows desktop path makes trouble otherwise.

% evalin('base', [   'print(gcf, ''-dpdf'', ''', [pathname, filename], ''');'    ]     );

evalin('base',  [ 'print(gcf,''-painters'',''-bestfit'',''-r300'',''-dpdf'' ,''' [ pathname filename ] ''');' ] ); %LM bugfix - save command was broken




if ~isempty(SD)
    SD.settings.lastDirEps = pathname;
    guidata(SD.figHandle, SD);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = optionsDialog(hObj,H,SD)
%%% Draw options dialog figure, wait for user 'OK' then set relevant fields in SD struct.  %%%
% Draw figure and panels %
hFig = figure('units', 'pixels','position',[(SD.screenPix(3)-500)/1.3 (SD.screenPix(4)-500)/1.3 500 500],'menubar','none','numbertitle','off','name','SCAn Options ..','closerequestfcn','delete(gcf)','renderer','painters');         
hPanAxesSize = uipanel('units','pixels','position',[15 320 240 150],'title','Figure Axes Size');
hPanDataLoad = uipanel('units','pixels','position',[15 220 240 80],'title','Data Loading');
uicontrol('units', 'pixels', 'position', [190 25 120 35], 'string', 'OK', 'callback', 'uiresume');
% Axes size panel %
panPropsTemp={'parent',hPanAxesSize,'units','normalized'};
lineTemp=0.7;   
uicontrol(panPropsTemp{1:end},'position',[0.05 lineTemp 0.5 0.2],'style','text','string','Fixed axes height (10ths screen height) ..');
uicontrol(panPropsTemp{1:end},'position',[0.6 lineTemp 0.25 0.2],'style','edit','tag','axesSizePlotHeight','BackgroundColor','w');
lineTemp=0.2; 
uicontrol(panPropsTemp{1:end},'position',[0.05 lineTemp 0.5 0.2],'style','text','string','Max Rows Per Figure ..');
uicontrol(panPropsTemp{1:end},'position',[0.6 lineTemp 0.25 0.2],'style','edit','tag','axesSizeMaxNRow','BackgroundColor','w');
% Data loading panel %
panPropsTemp={'parent',hPanDataLoad,'units','normalized'};
lineTemp=0.2;   
uicontrol(panPropsTemp{1:end},'position',[0.05 lineTemp 0.5 0.5],'style','text','string','Default cut tag ..');
uicontrol(panPropsTemp{1:end},'position',[0.6 lineTemp 0.25 0.5],'style','edit','tag','dataLoadDefaultCut','BackgroundColor','w');
%----%
uiwait
%----%
hUI = guihandles(hFig);
% Figure Axes Size % 
temp=get(hUI.axesSizePlotHeight,'string');  if ~isempty(temp);  SD.settings.fixedPlotHeight=str2double(temp);  end
temp=get(hUI.axesSizeMaxNRow,'string');     if ~isempty(temp);  SD.settings.maxRowNumber=str2double(temp);  end
% Data Load %
SD.settings.defaultCutTag = get(hUI.dataLoadDefaultCut,'string');
% Finish - assign and close %
guidata(SD.figHandle, SD);
close(hFig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = userVar(hObj,H,SD)
% UI for user to define trial or cell variable conditions. ('user' field of
% data structure).
% This function draws the GUI for defining user variables. Another subfunction 
% of SCAN_BASE, userVarRedraw, is responsible for 1) updating the UI to
% user actions, 2) setting the varaiables in the data struct.
% Note that userVarRedraw is called directly through a handle to the
% subfunction, rather than through the SCAN_BASE switchyard.
%%% Get relevant strings for figure %
level = get(hObj,'tag');   level = level(7:end);
switch level
case 'Data'
    item_list = SD.selData;
    var_list = {};
    for ii=1:length(SD.selData)
        data = evalin('base', SD.selData{ii});
        var_list = [var_list; fieldnames(data.user)];
    end
case 'Trial'
    data = evalin('base', SD.selData{1});
    item_list = getdatanames(data, 'trial', 'long', SD.selTrial);
    var_list = fieldnames(data.trials(1).user);               % All trials in a dataset should have same set of user variables set.
case 'Cell'
    data = evalin('base', SD.selData{1});
    item_list = getdatanames(data, 'cell', 'long',SD.selCell);
    var_list = fieldnames(data.trials(1).cells(1).user);      % All cells in a dataset should have same set of user variables set.
    dataUI.tet_list = cat(1, data.trials(1).cells.tet);
    dataUI.tet_list = dataUI.tet_list(SD.selCell);
end
var_list = unique(var_list);
var_list = [var_list; {''}];
%%% Draw figure %%%
screen = SD.screenChar(3:4);
lines_per_item = 2.5;
w_total = 80;
h_total = (length(item_list)*lines_per_item)+15;
size = [w_total h_total];
hFig = figure('name', 'SCAn: Define variable conditions', 'units', 'characters', 'position', [(screen - size)./2, size], ...
            'resize', 'on', 'defaultuicontrolunits', 'characters', 'menubar', 'none', 'numbertitle', 'off','windowstyle','modal');
% 'Select variable' section %
uicontrol('style', 'frame', 'position', [2 h_total-7.5 w_total-4 6]);
lbox = 8;
rbox = 45;
line1_y = h_total - 4;
line2_y = line1_y - 2;
uicontrol('style', 'text', 'position', [lbox line1_y 20 1.5], 'string', 'Choose variable .. ', 'horizontalalignment', 'left');
uicontrol('style', 'text', 'position', [rbox line1_y 20 1.5], 'string', '.. or define new.', 'horizontalalignment', 'left');
uicontrol('style', 'popupmenu', 'position', [lbox line2_y 30 1.5], 'string', var_list, 'value', length(var_list), ...
          'backgroundcolor', [1 1 1], 'tag', 'vrSel', 'callback', @userVarRedraw);
uicontrol('style', 'edit', 'position', [rbox line2_y 25 1.5], 'backgroundcolor', [1 1 1], 'callback', @userVarRedraw, 'tag', 'vrNew');
% Values %
for ii=1:length(item_list)
    line_y = (h_total-10.5)-((ii-1)*lines_per_item);
    uicontrol('style', 'text', 'position', [15 line_y 20 1], 'string', item_list{ii});
    % Value Edit boxes: convenient to have these handles in a vector array %  
    if strcmp(level, 'Cell')
        % In cell mode, Value boxes need callbacks for the 'lock_tet' feature %
        dataUI.values(ii) = uicontrol('style', 'edit', 'position', [45 line_y 20 1.5], 'backgroundcolor', [1 1 1], ...
                                      'callback', @userVarRedraw, 'tag', ['value' num2str(ii)]);
    else
        % Otherwise no callback %
        dataUI.values(ii) = uicontrol('style', 'edit', 'position', [45 line_y 20 1.5], 'backgroundcolor', [1 1 1]);
    end
end
% OK button etc. %
line_y = 2;
button_names = {'Done', 'Apply', 'Close'};
for ii=1:3
    uicontrol('style', 'pushbutton', 'position', [8+((ii-1)*23) line_y, 17 2], 'string', button_names{ii}, ...
              'callback', @userVarRedraw, 'tag', button_names{ii});
end
dataUI.lock_tet = 0;
dataUI.level = level;
guidata(hFig, dataUI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = userVarRedraw(hObj,eventData)
% Redraw function for userVar GUI. Importantly, also sets the values in the dataset structures.
% Note the input argument format - this function is called directly by
% 'userVar' subfunction callbacks (via f handles), not through main SCAN_BASE switchyard.
hUI = guihandles(hObj);
dataUI = guidata(hObj);
SD = gss;
%%% Get the variable name %%%
% To avoid ambiguities when setting .user values, will clear edit box whenever popupmenu is selected and visa versa %
switch get(hObj,'tag')
    case 'vrNew'
        var_set = get(hUI.vrNew, 'string');
        set(hUI.vrSel, 'value', length(get(hUI.vrSel, 'string')));    % Set 'value' to end (blank string)
    case 'vrSel'
        temp1 = get(hUI.vrSel, 'string');
        temp2 = get(hUI.vrSel, 'value');
        var_set = temp1{temp2};
        set(hUI.vrNew, 'string', '');
    case {'Done', 'Apply'}
        if ~isempty(get(hUI.vrNew, 'string'))
            var_set = get(hUI.vrNew, 'string');
        else
            temp1 = get(hUI.vrSel, 'string');
            temp2 = get(hUI.vrSel, 'value');
            var_set = temp1{temp2};
        end
end
%%% Either show values for selected variables, or set values in data structure %%%
tagObj = get(hObj,'tag');
switch strtok(tagObj, double('1234567890'))  % To protect against 'value12' type tags.
    case {'vrSel', 'vrNew'}
        % This will display any currently set values for the variable entered in the correct 'value' edit box %
        if isempty(var_set)
            return
        end
        switch dataUI.level
        case 'Data'
            for ii=1:length(dataUI.values)
                data = evalin('base', SD.selData{ii});
                try    val_set = data.user.(var_set);   catch;   val_set = '';   end
                set(dataUI.values(ii), 'string', val_set);
            end
        case {'Trial','Cell'}
            data = evalin('base', SD.selData{1});
            for ii=1:length(dataUI.values)
                if strcmp(dataUI.level, 'Trial')
                    u = data.trials( SD.selTrial(ii) ).user;
                else
                    u = data.trials(1).cells( SD.selCell(ii) ).user;
                end
                try    val_set = u.(var_set);   catch;   val_set = '';   end
                set(dataUI.values(ii), 'string', val_set);
            end
        end
    case {'Done', 'Apply'}
        % Look for the index of the variable entered, take the values entered, set in the structure %
        if isempty(var_set)
            errordlg('You need to define a variable (in one of the top boxes).');
            return
        end
        switch dataUI.level
        case 'Data'
            for ii=1:length(dataUI.values)
                data = evalin('base', SD.selData{ii});
                data.user(1).(var_set) = get(dataUI.values(ii),'string');
                assignin('base', SD.selData{ii}, data);
            end
        case {'Trial','Cell'}
            data = evalin('base', SD.selData{1});
            for ii=1:length(dataUI.values)
                if strcmp(dataUI.level, 'Trial')
                    data.trials( SD.selTrial(ii) ).user(1).(var_set) = get(dataUI.values(ii),'string');
                else
                    for jj=1:length(data.trials)    % For cells, set variable in every trial.
                        data.trials(jj).cells( SD.selCell(ii) ).user(1).(var_set) = get(dataUI.values(ii),'string');
                    end
                end
            end
            assignin('base', SD.selData{1}, data);
        end
        % Reset figure - include new var in popupmenu %
        s = get(hUI.vrSel, 'string');
        if ~any(strcmp(var_set, s))
            s{end} = var_set;
            s{end+1} = '';
        end
        set(hUI.vrSel, 'string', s);
        % Reset figure - set all to blank %
        set(dataUI.values, 'string', '');
        set(hUI.vrNew, 'string', '');
        set(hUI.vrSel, 'value', length(get(hUI.vrSel, 'string')));
        % If 'Done ' rather than 'Apply', close UI %
        if strcmp(tagObj, 'Done')
            delete(gcf);
        end
    case 'value'
        %%% This is the callback for the value edit boxes. Used to set values for other cells on tet %%%
        if dataUI.lock_tet
            val_ind = str2double(tagObj(6:end));
            v = get(dataUI.values(val_ind), 'string');
            set(dataUI.values( dataUI.tet_list==dataUI.tet_list(val_ind) ), 'string', v);
        end
    case 'Close'
        % Close figure, do nothing to data %
        delete(gcf);   
end









