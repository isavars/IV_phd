function [varargout] = load_set(fileName)
% Read inforation from a DACQ .set file.
%
%
%              [rtn] = load_set(trial_name);
%   [rtn setFileTxt] = read_set(trial_name);
%
%   rtn.duration
%   rtn.ADC_fullscale_mv
%   rtn.tracked_spots
%   rtn.lightBearing    e.g. [270 90 0 0]
%   rtn.colactive       e.g. [1 1 0 0]
%   rtn.gains(tet,ch)
%   rtn.eeg(x).channel  (x=1:nEEG)
%   rtn.eeg(x).gain
%   rtn.eeg(x).filter   (plus further filter spec for DacqUSB)
%
% setFileTxt is the unprocessed set file cell array, can be used with mTint
% KEY_VALUE.

% Check file exists, give a meaningful error if not %
fid=fopen([fileName '.set']);
if fid==-1
%     error('scan:load_set:setFileNotFound', ['Could not open set file ' fileName '.set']);
    ME = MException('scan:load_set:setFileNotFound', ['Could not open set file ' fileName '.set']);
    throw(ME);
else
    fclose(fid);
end

% Read file %
[key value] = textread([fileName '.set'], '%s %[^\n]');
txt = [cat(1,key) cat(1,value)];

% fid = fopen([fileName '.set']);
% txt = textscan(fid, '%s');
% fclose(fid);

%%% These are straightforward one value fields %%%
rtn.tracked_spots = str2double(getValue(txt, 'tracked_spots'));
rtn.xmin = str2double(getValue(txt, 'xmin'));
rtn.xmax = str2double(getValue(txt, 'xmax'));
rtn.ymin = str2double(getValue(txt, 'ymin'));
rtn.ymax = str2double(getValue(txt, 'ymax'));
rtn.sw_version = getValue(txt, 'sw_version');   % Don't STR2DOUBLE - can be 4.00a
rtn.trial_time = getValue(txt, 'trial_time');
rtn.ADC_fullscale_mv = str2double(getValue(txt, 'ADC_fullscale_mv'));
temp = getValue(txt, 'duration'); rtn.duration = str2double(temp(1:end));%IV removed -1 causing problems for trials with more than 3 digit durations
% Identify the generation of DACQ %
if rtn.ADC_fullscale_mv == 1500
    rtn.dacq_version = 'usb';
elseif strcmp(rtn.sw_version(1),'1')
    rtn.dacq_version = 'legacy_1';
elseif strcmp(rtn.sw_version(1),'2') || strcmp(rtn.sw_version(1),'3')
    rtn.dacq_version = 'legacy_2';  % DACQ2 sw_version numbers start '2' and '3'
end

% Trial Duration %
if strcmp(rtn.dacq_version, 'legacy_1')
    temp = getValue(txt, 'duration'); % DACQ1 stores '900 s'
    rtn.duration = str2double(temp(1:end));%IV removed -1 causing problems for trials with more than 3 digit durations
else
    rtn.duration = str2double(getValue(txt, 'duration'));
end

% light parameters: make a 1:4 vector %
lightParams = {'lightBearing', 'colactive'};
for ii=1:2
    for jj=1:4
        rtn.(lightParams{ii})(jj) = str2double(getValue(txt, [lightParams{ii} '_' num2str(jj)]));
    end
end

% Gains %
% for ii=1:16
%LM bugfix
for ii=1:32
    for jj=1:4
        rtn.gains(ii,jj) = str2double(getValue(txt, ['gain_ch_' num2str((ii-1)*4 + jj-1)]));
    end
end
rtn.fullscale = (rtn.ADC_fullscale_mv ./ rtn.gains) .* 1000;

%%% EEG Channels %%%
% Which channels EEGs recorded? %
% ch = zeros([1 16]);
recordingChannel = zeros([1 64]); %LM bugfix
if strcmp(rtn.dacq_version,'usb')
    % DacqUSB %
%     if ~str2double(getValue(txt,'nullEEG'))
%         for ii=1:16 %% There are occasional DACQs which only have EEG 1-4. This loop therfore needs to fail gracefully.
        %LM bugfix
        for ii=1:length(recordingChannel) %% There are occasional DACQs which only have EEG 1-4. This loop therfore needs to fail gracefully.
            temp=getValue(txt,['saveEEG_ch_' num2str(ii)]);
            if isempty(temp)
                break;  
            elseif str2double(temp) % temp is '1' or '0' for EEG used or not.
                recordingChannel(ii) = str2double(getValue(txt,['EEG_ch_'  num2str(ii)]));
            end
        end
%     end
elseif strcmp(rtn.dacq_version,'legacy_2')
    % Dacq 2 %
    if ~str2double(getValue(txt,'nullEEG'))
        for ii=1:2
            if str2double(getValue(txt,['saveEEG_ch_' num2str(ii)]))
                recordingChannel(ii) = str2double(getValue(txt,['EEG_ch_'  num2str(ii)]));
            end
        end
    end
elseif strcmp(rtn.dacq_version,'legacy_1')
    % Dacq1 %
    for ii=1:8
        if str2double(getValue(txt, ['EEGmap_' num2str(ii)]));
            recordingChannel(ii) = ii*4;
        end
    end
end
EEGSlotActive          = find(recordingChannel);                  % Now a list of eeg channels in use ..
recordingChannel = recordingChannel(recordingChannel~=0);   %  .. and their recording slots. (recordingChannel=0 if slot not in use).
% Get signal source channel, gain, filters %
if isempty(recordingChannel)
    rtn.eeg = []; % In case of null EEG
else
    for ii=1:length(recordingChannel)
        mode = str2double(getValue(txt, ['mode_ch_' num2str(recordingChannel(ii)-1)]));
        if any(mode==[1 3]) % Mode B or -B
            %%LM: update (30/04/2015) - when using the new free referencing system the
            %%old code doesn't load the eeg correctly ('modeanalog32' == 1 is free
            %%referencing)
            if str2double(getValue(txt,'modeanalog32')) == 1
                chTemp = getValue(txt, ['b_in_ch_' num2str(recordingChannel(ii)-1)]);
            else
                chTemp = getValue(txt, ['b_in_ch_' num2str(recordingChannel(ii)-1)]);
                if strcmp(rtn.dacq_version,'usb') % DacqUSB-there is a secondary reference to the actual ch.
                    chTemp = getValue(txt, ['ref_' chTemp]);
                end
            end
            rtn.eeg(ii).channel = str2double(chTemp) + 1;
        else
            rtn.eeg(ii).channel = recordingChannel(ii);
        end
        rtn.eeg(ii).recordingChannel = recordingChannel(ii);       % Actual recording channel, not signal source (useful when approaching raw data).
        rtn.eeg(ii).EEGSlot          = EEGSlotActive(ii);  % Is this EEG1, EEG2, etc (according to the 'record>setup>EEG' tab.
        rtn.eeg(ii).scalemax         = rtn.fullscale( ceil(recordingChannel(ii)/4), recordingChannel(ii)-(((ceil(recordingChannel(ii)/4))-1)*4) );
        rtn.eeg(ii).filter           = str2double(getValue(txt, ['filter_ch_' num2str(recordingChannel(ii)-1)]));
        if strcmp(rtn.dacq_version,'usb')
            % DacqUSB only - fuller filter spec %
            fSpec = {'filtresp', 'filtkind', 'filtfreq1', 'filtfreq2', 'filtripple'};
            for jj=1:length(fSpec)
                rtn.eeg(ii).(fSpec{jj}) = str2double(getValue(txt, [fSpec{jj} '_ch_' num2str(recordingChannel(ii)-1)]));
            end
        end
    end
end

% Output %
varargout{1} = rtn;
if nargout==2
    varargout{2} = txt;
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rtn] = getValue(txt,keyStr)
% Get a value from value-key cell array %
ind = strcmp(keyStr,txt(:,1));
if sum(ind) == 0
    rtn = [];
else
    rtn = txt{ind,2};
end
