 
function [varargout] = load_main(SD, varargin)

% Load a dataset, consisting of data over a series of trials.
%
%       (1) load_main(SD);
%       (2) data = load_main(SD, LS, setFileNames, pathName);
%
% In mode (1) all parameters are gathered via GUI, No need for return
% argument - will assign loaded data to workspace.
%
% Mode (2) allows you to call LOAD_MAIN functionally. (Will one day be 
% utilized by some sort of 'reload' function). For 'LS' argument, see 
% LOAD_SELECTDATA.
%
% 'setFileNames' - cell array of strings with names of set files for trials to load.
% 'pathName' - path to directory with data.

dirReset = pwd;
if isempty(varargin)
    %%% Get parameters via GUI %%%
    % Go to last data directory %
    try
        cd(SD.settings.lastDirLoad);
    catch
        SD.settings.lastDirLoad = pwd;
    end
    %%% Get Trial names %%%
    [trialNamesTemp, pathName] = uigetfile({'*.set'}, 'Select Trials to Load .. ', 'multiselect', 'on');
    if ~ischar(pathName);   return;   end
    if ischar(trialNamesTemp);   
        trialNames{1} = trialNamesTemp;     % UIGETFILE output is char for single selections ..
    else
        trialNames = trialNamesTemp;        % .. cell string otherwise.
    end
    SD.settings.lastDirLoad = pathName;
    trialNames = sort(trialNames);
    %%% Get Cell Selection and Data Parameters %%%
    LS = load_selectdata(trialNames, pathName);  % LS = Loading Selection
else
    %%% Use parameters supplied by calling function %%%
    LS = varargin{1};
    trialNames = varargin{2};
    pathName = varargin{3};
end

% Initialise data struct %
data = structInit;
data.trials(2:length(trialNames)) = data.trials(1); % So all trials have default .user fields.

LS.pathName = pathName;       % Save these to data struct for convenience
LS.trialNames = trialNames;   % if reloading data.
data.load_selection = LS;     %

%%% Load Data %%%
cd(pathName);
for ii=1:length(trialNames)
    if LS.verboseLoading
        fprintf(1,'\nLoading Data : %s; ', trialNames{ii}(1:end-4));
    end
    
    %%% Set file %%%
    [setFile setFileText] = load_set(trialNames{ii}(1:end-4));
    
    % Truncate overhanging 1s. %
    rv = 10;  % Clip anything not a multiple of 10s.
    if LS.clipOverhang && rem(setFile.duration, rv)~=0
        disp(['!!NOTE!! Overhanging 1s? Trial duration is ' num2str(setFile.duration) ...
              ' - truncating to ' num2str(setFile.duration - rem(setFile.duration, rv)) '.']);
        setFile.duration = setFile.duration - rem(setFile.duration, rv);
        setFileText{strmatch('duration', setFileText(:,1)),2} = num2str(setFile.duration);  % For mTint functions.
    end

    %%% Positions %%%
    pos = load_pos(trialNames{ii}(1:end-4), setFile.tracked_spots);  % Always need this, as it gives the most reliable trial duration.
    if LS.loadPos
        if LS.posMode==1
            % Load from .dat file %
            dat = load_dat(trialNames{ii}(1:end-4));
            data.trials(ii).x = int16(dat.x + 1);   % To shift origin pixel from
            data.trials(ii).y = int16(dat.y + 1);   % 0,0 to 1,1
            data.trials(ii).dir = int16(dat.dir);
            data.trials(ii).speed = single(dat.speed);
            data.trials(ii).ppm = dat.ppm;
            data.trials(ii).sample_rate = dat.sample_rate;
            data.trials(ii).window_x = 512;
            data.trials(ii).window_y = 512;
        elseif LS.posMode==2
            %%% Load from DACQ files %%%
            % If not enough pos samples, reduce the stated trial duration (DACQ1 + esc'd trials).
            if size(pos.led_pos,1) < setFile.duration*pos.sample_rate
                setFile.duration = floor(size(pos.led_pos,1)/pos.sample_rate);
                setFileText{strmatch('duration', setFileText(:,1)),2} = num2str(setFile.duration);  % For mTint functions.
            end
            % Clip overhanging 1s, or other extraneous pos samples. (FLOOR because DACQ1 SR = 46.875Hz).
            pos.led_pos = pos.led_pos(1:(floor(pos.sample_rate*setFile.duration)),:,:);    
            pos.led_pix = pos.led_pix(1:(floor(pos.sample_rate*setFile.duration)),:,:);
            pos.header{strmatch('duration', pos.header(:,1)),2} = num2str(setFile.duration);      % For mTint functions.
            % Scale, if required %
            if ~isnan(LS.posScale)
                pos.led_pos = floor(pos.led_pos .* (LS.posScale/pos.pixels_per_metre));
                data.trials(ii).original_ppm = pos.pixels_per_metre;
                pos.pixels_per_metre  = LS.posScale;
                pos.header{strmatch('pixels_per_metre',pos.header(:,1)),2} = num2str(LS.posScale);  % For POSTPROCESS_POS_DATA
            end
            % Interpolate, Filter and Smooth %
            [xy, dir, speed] = postprocess_pos_data(pos, LS.posMaxSpeed, LS.posSmooth/1000, setFileText, LS.posHead);  % Smooth window: convert ms to s.
            data.trials(ii).x = int16( floor(xy(:,1)) +1 );
            data.trials(ii).y = int16( floor(xy(:,2)) +1 );
            if size(dir,1)<size(dir,2);  dir=dir';  end  % POSTPROCESS_POS_DATA gives dir in different orientation for one or two lights.
            data.trials(ii).dir = single( round(dir.*10)./10 );
            data.trials(ii).speed = single(speed);
            data.trials(ii).sample_rate = pos.sample_rate;
            data.trials(ii).ppm = pos.pixels_per_metre;
            data.trials(ii).window_x = setFile.xmax - setFile.xmin; % Can be inaccurate in pos file header.
            data.trials(ii).window_y = setFile.ymax - setFile.ymin; % Also note no scaling: like TINT system.
        end
    end
    
    %%% Spike Times, Waveforms %%%
    if LS.loadSpikes
        firstTet = 1;
        for jj=1:length(LS.tets)
            if ~isempty(LS.tets(jj).cellList)
                if LS.verboseLoading;    fprintf(1,'Tet %d, ',jj);    end
                cellsTemp = load_tet(trialNames{ii}(1:end-4), LS, jj, setFile, pos);  % 'pos' passed so as to know 'accurate' trial dur for DACQ1 
                if firstTet
                    data.trials(ii).cells(1:length(cellsTemp)) = cellsTemp;
                else
                    data.trials(ii).cells(end+1:end+length(cellsTemp)) = cellsTemp;
                end
                firstTet = 0;
            end
        end
    end
    
    %%% EEG %%%
    if any([LS.loadEEG1 LS.loadEEGAll LS.loadEGF])
        if isempty(setFile.eeg)
            data.trials(ii).eeg = [];  % In case of null EEG
        else
            for jj=1:length(setFile.eeg)
                % What's the filename? - Depends on Dacq version %
                if strcmp(setFile.dacq_version, 'usb')
                    % DacqUSB %
                    if jj==1;   ext = 'eeg';   else   ext = ['eeg' num2str( setFile.eeg(jj).EEGSlot )];   end
                else
                    % Dacq 1, 2 %
                    if jj==1;   ext = 'eeg';   else   ext = 'eg2';   end
                end
                eegTemp = load_eeg([trialNames{ii}(1:end-3) ext]);
                data.trials(ii).eeg(jj).eeg = eegTemp.eeg;
                % Below: nPos*5 is more reliable than EEG_sampRt for DACQ1 sampRt (46.875Hz) + nasty trial durations, if pos data exists
                if ~isempty(data.trials(ii).x) && ~isempty(data.trials(ii).eeg(jj).eeg)
                    data.trials(ii).eeg(jj).eeg = data.trials(ii).eeg(jj).eeg(1:length(data.trials(ii).x)*5);
                end
                data.trials(ii).eeg(jj).sample_rate = eegTemp.sample_rate;
                % EEG metadata from setFile %
                eegData = fieldnames(setFile.eeg);
                for kk=1:length(eegData)
                    data.trials(ii).eeg(jj).(eegData{kk}) = setFile.eeg(jj).(eegData{kk});
                end
                if ~LS.loadEEGAll %
                    break         % Load only EEG1
                end               %
            end
        end
    end
    
    data.trials(ii).trialname = trialNames{ii}(1:end-4);
    data.trials(ii).rat = pathName;
    data.trials(ii).dur = setFile.duration;

    %%% Log file, Stm file etc %%%
    data.trials(ii).stim_times=load_stm(trialNames{ii}(1:end-4));
    fid = fopen([trialNames{ii}(1:end-4), '.log']);
    if fid ~= -1
        temp=textscan(fid,'%s');
        data.trials(ii).log = temp{1};
        fclose(fid);       
    end  
end
cd(dirReset);

% Get rid of cell with no spikes in any trial (implements 'all with spikes' option) %
data = trimNullCells(data, LS);

if nargout == 0
    % Assign output in workspace %
    varargout = {};
    assignin('base', LS.dataName, data);
    guidata(SD.figHandle, SD);
    scan_base('listRefreshData',SD.figHandle);
else
    % Assign output to calling function %
    varargout{1} = data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rtn] = structInit()
%%% Initialise dataset structure %%%
% Whole structure tags %
rtn.is_scan = 1;
rtn.rate_maps = struct([]);
rtn.user.day = '';
rtn.user.age = '';
rtn.user.brain_region = '';
rtn.user.genotype = '';
% General trial info tags %
rtn.trials.rat = [];
rtn.trials.trialname = [];
rtn.trials.ppm = [];
rtn.trials.original_ppm = [];
rtn.trials.dur = [];
rtn.trials.sample_rate = [];
rtn.trials.log = [];
rtn.trials.user.environment = '';
rtn.trials.user.shape = '';
rtn.trials.user.n_exp = '';
% Trial Position Data %
rtn.trials.x = [];
rtn.trials.y = [];
rtn.trials.speed = [];
rtn.trials.dir = [];
% Trial EEG Data %
rtn.trials.eeg(1).eeg = [];
rtn.trials.eeg(1).sample_rate = [];
% Cells: Info % 
%%% IMPORTANT: field names and declaration order should match LOAD_TET %%%
rtn.trials.cells(1).st = [];
rtn.trials.cells(1).wf_all = [];
rtn.trials.cells(1).wf_means = [];
rtn.trials.cells(1).wf_amps = [];
rtn.trials.cells(1).cut = [];
rtn.trials.cells(1).cellnum = [];
rtn.trials.cells(1).tet = [];
rtn.trials.cells(1).scalemax = [];
rtn.trials.cells(1).wf_mode = 0;
rtn.trials.cells(1).user.brain_region = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = trimNullCells(data, LS)
% Remove from dataset cells with no spikes in any trial. (Implements 'all
% with spikes' option - cells 1-19 loaded, then spikeless cells removed).
% If no spike data selected %
if ~LS.loadSpikes;   return;   end
% Get list of relevant tets %
trimList = cat(1, LS.tets.trimNull);
if sum(trimList)==0;   return;   end  
% Get a list of cells which are on the tets with trim_null=1 %
trimCellInd = [];
for ii=1:length(data.trials(1).cells)
    if trimList( data.trials(1).cells(ii).tet )  % Tet num is index into trimList (logical)
        trimCellInd = [trimCellInd; ii];
    end
end
% Out of trimCellInd, now look for cells null across all trials %
nullCellsInd = [];
for ii=1:length(trimCellInd)
    for jj=1:length(data.trials)
        if length(data.trials(jj).cells(trimCellInd(ii)).st) > 0
            break
        end
        if jj == length(data.trials)                        % if all trials for this cell are no spikes,
            nullCellsInd = [nullCellsInd, trimCellInd(ii)]; % mark as null.
        end                                                 %
    end
end
% Trim the data %
notNullInd = setdiff(1:length(data.trials(1).cells), nullCellsInd);    % this is the cells we want to keep
for ii=1:length(data.trials)
    data.trials(ii).cells = data.trials(ii).cells(notNullInd);
end 
















