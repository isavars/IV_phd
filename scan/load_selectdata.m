function [Rtn] = load_selectdata(trialNames, pathName)
% UI for data load selection.
%
%        [Rtn] = load_selectdata(trialNames, pathName);
%
% Rtn.dataName
% Rtn.loadPos
% Rtn.loadSpikes
% Rtn.loadEEG1
% Rtn.loadEEGAll
% Rtn.loadEGF
% Rtn.posMode
% Rtn.wfMode
% Rtn.posMaxSpeed
% Rtn.posSmooth
% Rtn.posHead
% Rtn.posScale
% Rtn.clipOverhang
% Rtn.verboseLoading
% Rtn.skipMissingTets
% 
% Rtn.tets(1:8).cellList
% Rtn.tets(1:8).trimNull
% Rtn.tets(1:8).cutTag
% Rtn.tets(1:8).cut
% Rtn.tets(1:8).load0

%%% Generate default namer for data struct %%%
animalName = fliplr(strtok(fliplr(pathName), '/\'));            % Hoping lowest folder of pathname is animal name
defaultDataName = [animalName, '_', trialNames{1}(1:end-5)];    % Hoping to trim off trial letter (+ .set)
% Is first letter number? - put 'r' at the front %
if ~isempty(strfind('0123456789', defaultDataName(1)))
    defaultDataName = ['r', defaultDataName];
end
% Check not overwriting other data %
existingData = evalin('base', 'who');
n=1;
while ~isempty(strmatch(defaultDataName, existingData))
    defaultDataName = [defaultDataName, '_', num2str(n)];
    n=n+1;
end

Rtn = [];
SD.screenPix = [1 1 1280 800];
%%% Draw figure %%%
hFig = figure('units', 'pixels', 'position', [(SD.screenPix(3)-1065)/1.3 (SD.screenPix(4)-510)/1.7 1065 510], ...
               'menubar','none','numbertitle','off','name','Select data to load ..','closerequestfcn','delete(gcf)','renderer','painters');
% Panels %
hPanName = uipanel('units','pixels','position',[5 415 290 90],'title','Data Name');
hPanData = uipanel('units','pixels','position',[5 5 290 405],'title','Data Loading Options');
hPanCells = uipanel('units','pixels','position',[305 75 750 430],'title','Select Cells','tag','selectCells');
% OK, Cancel Buttons %
uicontrol('units', 'pixels', 'position', [550 20 120 35], 'string', 'OK', 'callback', 'uiresume');
uicontrol('units', 'pixels', 'position', [700 20 120 35], 'string', 'Cancel', 'callback', 'uiresume');

%%% Data Name Panel %%%
uicontrol('parent', hPanName, 'units', 'normalized', 'position', [0.1 0.6 0.8 0.3], 'style', 'text', 'string', 'Enter name for this data set:', 'fontsize', 10, 'fontweight', 'demi', 'horizontalalignment', 'left');
uicontrol('parent', hPanName, 'units', 'normalized', 'position', [0.1 0.3 0.8 0.3], 'style', 'edit', 'backgroundcolor', 'w', 'string', defaultDataName, 'horizontalalignment', 'left','tag','dataName');

%%% Data Options Panel %%%
ySpacer = [0.03:0.06:0.15 0.24:0.07:0.94];
indent = [0.05 0.1 0.15];
% Main options list %
     %'labelString','UIStyle',indent,'tag',defaultValue; ...
cv = {'Load EGFs','checkbox',1,'egfCheck',0; ...
      'All recorded EEGs','checkbox',1,'eegAllCheck',1; ...
      'EEG1 only','checkbox',1,'eeg1Check',1; ...
      'All Waveform Samples','radiobutton',2,'wfMode2',0;...
      'Amps + Mean Waveforms','radiobutton',2,'wfMode1',1;...
      'Spike Times Only','radiobutton',2,'wfMode0',0;...
      'Spikes','checkbox',1,'spkCheck',1;...
      'Max Speed (m/s)','text',3,'posMaxSpeed',0;...
      'Smooth Window (ms)','text',3,'posSmooth',0;...
      'Scale (to this ppm):','text',3,'posScale',0;...
      'Head position:','text',3,'posHead',0;...
      'Load from DACQ files','radiobutton',2,'posMode2',1;...
      'Load from .dat file','radiobutton',2,'posMode1',0;...
      'Positions','checkbox',1,'posCheck',1};
for ii=1:length(ySpacer)
    hTemp = uicontrol('parent',hPanData,'units','normalized','position',[indent(cv{ii,3}) ySpacer(ii) 0.8 0.04],'value',cv{ii,5},...
        'style', cv{ii,2}, 'string', cv{ii,1},'horizontalalignment','left','callback',@fRedraw,'tag',cv{ii,4});
    if ii==1                         % Note, 05/09/08. EGF loading checkbox
        set(hTemp, 'enable', 'off');  % disabled until EGF loading is sorted out.
    end                               
end
% DACQ position loading edit boxes %
defaultValue = {'4','400', '', '0.5'};   % Max Speed (m/s), Smooth Window (ms), Scale (to this ppm), Head position.
yPosInd = 8:11;
for ii=1:length(yPosInd)
    hTemp = uicontrol('parent', hPanData, 'units', 'normalized', 'position', [0.6 ySpacer(yPosInd(ii)) 0.2 0.04], ...
              'style', 'edit', 'string', defaultValue{ii},'tag',cv{yPosInd(ii),4},'backgroundcolor','w');       
end

%%% Cells Panel %%%
rowPos = fliplr(0.05:0.05:0.82);  % Flip to go top to bottom
titleRowPos = 0.92;
colPos = [0.015 0.07 0.155 (0:0.039:0.16)+0.2  (0:0.039:0.16)+0.4 (0:0.039:0.16)+0.6 (0:0.039:0.12)+0.8 0.96];  
% Columns for each tet: label, AWS box, regular cell boxes and cell box frames %
for ii=1:16
    uicontrol('parent', hPanCells, 'style', 'text', 'string', ['Tet ' num2str(ii)], 'units', 'normalized', 'position', [colPos(1) rowPos(ii) 0.06 0.04],'horizontalalignment','left');
    cellTicks(ii,1) = uicontrol('parent',hPanCells,'style','checkbox','units','normalized','position', [colPos(3) rowPos(ii) 0.025 0.06],'callback',@fRedraw,'tag',['allWithSpikes' num2str(ii)]);
    for jj=2:21
        cellTicks(ii,jj) = uicontrol('parent',hPanCells,'style','checkbox','units','normalized','position',[colPos(jj+2) rowPos(ii) 0.025 0.06]);
    end
end
% Cut tag boxes done all together, so they all together in tab order of UI.
for ii=1:16
    cutTag(ii) = uicontrol('parent', hPanCells,'style', 'edit', 'string', ['_' num2str(ii)], 'units', 'normalized', 'position', [colPos(2) rowPos(ii) 0.06 0.05], 'backgroundcolor', [1 1 1]);
end
set(hPanCells,'userdata',cellTicks);  % So fRedraw function can access nice handle array
% Column Labels %    
uicontrol('parent', hPanCells,'style', 'text', 'string', {'Cut', 'Indentifier'} , 'units', 'normalized', 'position', [colPos(2)-0.005 titleRowPos 0.07 0.06]);
uicontrol('parent', hPanCells, 'style', 'text', 'string', 'All With Spikes', 'units', 'normalized', 'position', [colPos(3)-0.01 titleRowPos-0.03 0.045 0.1]);
for ii=1:20
    uicontrol('parent', hPanCells,'style', 'text', 'string', num2str(rem(ii,20)), 'units', 'normalized', 'position', [colPos(ii+3) titleRowPos 0.02 0.03]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uiwait;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If 'Cancel', Rtn is [] %
if strcmp(get(gco,'string'),'Cancel')
    Rtn = [];   delete(hFig);   return
end

%%% Wrap up selections into output structure and return %%%
% Cells Selection for each Tet %
% NOTE: this code also checks whether all required cut files exist - use subfunction in a WHILE loop so that figure
% is held open until missing cut files are sorted out - without this, missing cuts are a big pain in the arse.
[Tets missingList] = getCellValues(cellTicks, cutTag, trialNames, pathName);
while ~isempty(missingList)
    fprintf(1,'\n SCAN Load Data: The following cut files could not be read:\n');   disp(missingList);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uiwait;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Tets missingList] = getCellValues(cellTicks, cutTag, trialNames, pathName);
end
Rtn.tets = Tets;
% Data Options %
hUI = guihandles(hFig);
Rtn.dataName = get(hUI.dataName, 'string');
Rtn.loadPos = get(hUI.posCheck, 'value');
Rtn.loadSpikes = get(hUI.spkCheck, 'value');
Rtn.loadEEG1 = get(hUI.eeg1Check, 'value');
Rtn.loadEEGAll = get(hUI.eegAllCheck, 'value');
Rtn.loadEGF = get(hUI.egfCheck, 'value');
Rtn.posMode = find(cell2mat(get([hUI.posMode1 hUI.posMode2],'value')));
Rtn.wfMode = find(cell2mat(get([hUI.wfMode0 hUI.wfMode2 hUI.wfMode1],'value'))) - 1; % Note: mode1 = all data, mode2 = wf+amps, for back compat.
Rtn.posMaxSpeed = str2double(get(hUI.posMaxSpeed(1),'string'));
Rtn.posSmooth = str2double(get(hUI.posSmooth(1),'string'));
Rtn.posHead = str2double(get(hUI.posHead(1),'string'));
Rtn.posScale = str2double(get(hUI.posScale(1),'string'));
Rtn.clipOverhang = 1;           % Need to write UI for this.
Rtn.verboseLoading = 1;         % No UI for this, assume that loading from GUI is always verbose (ie load_main reports trials and tet fiels read).
Rtn.skipMissingTets = 0;        % If set to 1, missing cut files will be 'skipped over' - no error thrown, cells that should have been specified are in structure, but with no spikes. A warning will be given. This also isn't in the UI, only used for automated load-ins.
delete(hFig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tets missingList] = getCellValues(cellTicks, cutTag, trialNames, pathName)
% Subfunction to 1) roll up 'cells' pane of UI into return structure, 2) check to see 
% whether all required cut files exist. This code is sectioned off as a subfunction so 
% that it can be called in a WHILE loop until user has sorted out all cut files.
for ii=1:16
    if get(cellTicks(ii,1),'value')
        Tets(ii).cellList = (1:30)';    % This is for 'All with spikes' option
        Tets(ii).trimNull = 1;          %
    else
        Tets(ii).cellList = find(cell2mat(get(cellTicks(ii,1:end-1),'value'))) - 1; % cellTicks(ii,1:end-1) - don't get cell 0 value.
        Tets(ii).trimNull = 0;
    end
    Tets(ii).load0 = get(cellTicks(ii,end),'value');
    Tets(ii).cutTag = get(cutTag(ii), 'string'); 
    Tets(ii).cut = '';  % This field is for LOAD_TET to check against - for RELOAD type functions, will contain a cut file name.
end
missingList = {};
for ii=1:16
    if ~isempty(Tets(ii).cellList)
        for jj=1:length(trialNames)
            fid = fopen([pathName trialNames{jj}(1:end-4) Tets(ii).cutTag '.cut'], 'r');
            if fid==-1
                missingList{end+1} = [trialNames{jj}(1:end-4) Tets(ii).cutTag '.cut'];
            else
                fclose(fid);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = fRedraw(hObj,eventData)
% This subfunction redraws the UI to take account of selection contingent
% enabling/disabling and radio button groups.
hUI = guihandles(hObj);
switch get(hObj,'tag')
    case 'posCheck'
        if get(hObj,'value')
            set([hUI.posMaxSpeed hUI.posSmooth hUI.posHead hUI.posScale hUI.posMode1 hUI.posMode2], 'enable', 'on');
        else
            set([hUI.posMaxSpeed hUI.posSmooth hUI.posHead hUI.posScale hUI.posMode1 hUI.posMode2], 'enable', 'off');
        end
    case 'posMode1'
        set([hUI.posMaxSpeed hUI.posSmooth hUI.posHead hUI.posScale], 'enable', 'off');
        set(hUI.posMode2, 'value', 0);
        set(hUI.posMode1, 'value', 1);
    case 'posMode2'
        set([hUI.posMaxSpeed hUI.posSmooth hUI.posHead hUI.posScale], 'enable', 'on');
        set(hUI.posMode1, 'value', 0);
        set(hUI.posMode2, 'value', 1);
    case 'spkCheck'
        if get(hObj,'value')
            set([hUI.wfMode0 hUI.wfMode1 hUI.wfMode2 get(hUI.selectCells,'children')'], 'enable', 'on');
        else
            set([hUI.wfMode0 hUI.wfMode1 hUI.wfMode2 get(hUI.selectCells,'children')'], 'enable', 'off');
        end
    case {'wfMode0','wfMode1','wfMode2'}
        set([hUI.wfMode0 hUI.wfMode1 hUI.wfMode2], 'value', 0);
        set(hObj, 'value', 1);
    case {'allWithSpikes1','allWithSpikes2','allWithSpikes3','allWithSpikes4','allWithSpikes5','allWithSpikes6','allWithSpikes7','allWithSpikes8'}
        cellTicks = get(get(hObj,'parent'),'userdata');
        str = get(hObj,'tag');   hCells = cellTicks(str2double(str(end)), 2:end);
        if get(hObj,'value')
            set(hCells,'enable','off');
        else
            set(hCells,'enable','on');
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Should put this in here, if I'm going to put it at all.

% % If loading from .dat files, Check if .dat files exist for all the trials %
% cd(data_dir);
% if (data_types.pos == 1) || (data_types.eeg == 1)
% 	dat_exist = zeros(1,length(trial_names));
% 	for i=1:length(trial_names);
%         for j=1:16  % check every tetrode
%             fid = fopen([trial_names{i}(1:end-4), '_', num2str(j),'.dat'], 'r');
%             if fid ~= -1;
%                 fclose(fid);
%                 dat_exist(i) = 1;
%                 break   
%             end
%         end
% 	end
% 	if any(dat_exist == 0)
%         no_dat = find(dat_exist == 0);
%         fprintf(1,'LOADDATA aborted, because no .dat file for trial(s): \n');
%         for i=1:length(no_dat)
%             fprintf(1,'%s\n', trial_names{no_dat(i)}(1:end-4));
%         end
%         cd(dir_reset);
%         return   
% 	end
% end








