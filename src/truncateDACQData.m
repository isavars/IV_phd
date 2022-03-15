function truncateDACQData( timeS,cutFileFlag,path2file )
% This function will truncate DACQ data between timeS(1) and time S(2). All
% trial data will be truncated. This is useful when e.g. at some point the 
% headstage unplugged and one wants to salvage the data recorded prior to 
% this or when only certain periods of the trial, based e.g. on animal
% behaviour should be extracted.
% Data is selected via UI or path is supplied by 'path2file'.
%
%
%  Usage:   truncateDACQData( timeS ) 
%           truncateDACQData( timeS,cutFileFlag )
%           truncateDACQData( timeS,cutFileFlag, path2file)
%           truncateDACQData( timeS,[], path2file )
%
%  Inputs:  
%           timeS       - [startTimes stopTimes] array or [stopTime] in seconds. In
%                         latter case we will assume startTime = 0
%           cutFileFlag - true/false (default) - indicates whther cutfiles
%                         should be truncated as well
%           path2file   - full path to file (including name & extension)
%         
%
%  LM 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% pre stuff
defaultDir = 'Z:\lmuessig\!postDoc\recording_data\CA1\'; % 'Z:\lmuessig\!postDoc\recording_data\CA1\'; 'D:\tempAdultData\';
if ~isdir(defaultDir)
    defaultDir = 'C:\';
end

%% input parseing
% sanity check
if any(diff(timeS,[],2)) < 0
    error('Time in every time chunk has to increase my friend. This code will terminate here!');
end
% if only stop time is supplied
if length(timeS) == 1 
    timeS = [0 timeS];
end
% default cut file flag is false
if nargin == 1 || isempty(cutFileFlag)
    cutFileFlag = false;
end

%% select data
% grab data
if nargin < 3
    [setName, dataDir] = uigetfile([defaultDir '*.set'],'Select .set File of Data To Truncate'); % UI
    if ~ischar(dataDir)
        warning('Loading cancelled. This code terminates here!');
        return
    end
else
    % get input in same format as if selected from UI 
    [dataDir,setName,ext] = fileparts(path2file);
    
    dataDir = [dataDir '\'];
    setName = [setName ext];
    clear ext
end

% backup dir - origina data will be moved here
backUpDir = [dataDir 'original_data\'];
if ~isdir(backUpDir)
    mkdir(backUpDir);
end
% dum = dir(backUpDir);
% if ~isempty(backUpDir)
%     
% end
% clear dum
% %% change set set file
[~,setNameNoExt,~] = fileparts([dataDir setName]); %remove extension

[~, sFileText] = load_set([dataDir setNameNoExt]); %load set files (use Scan function for convenience)
durInd = strcmp(sFileText(:,1),'duration');
trialDurOrg = sscanf(sFileText{durInd,2},'%d');
trialDurNew = sum(diff(timeS,[],2));

%sanity check
if trialDurNew > trialDurOrg
    error(['There is something wrong with the times you supplied. The truncated time (' num2str(trialDurNew) ') would be longer than the original' ...
                'trial time (' num2str(trialDurOrg) ')']) 
end
sFileText{durInd,2} = num2str(trialDurNew); %update trial duration
    
% change a few fields in header
timeInd = strcmp(sFileText(:,1),'trial_time');
trialTime = sFileText(timeInd,2);
trialTimeNew = datetime(trialTime{:},'inputformat','HH:mm:ss','format','HH:mm:ss') + seconds(timeS(1,1));
sFileText{timeInd,2} = char(trialTimeNew);
% back up original file
movefile([dataDir setName],[backUpDir setName]); 
% % write file
fid = fopen([dataDir setName],'w');
for j = 1:length(sFileText)
    fprintf(fid,'%s %s\r\n',sFileText{j,1},sFileText{j,2});
end
fclose(fid);
fprintf('Wrote .set file for truncated data: %s%s\n##########\n',dataDir, setName);
% clear timeInd durInd trialTime 

%% truncate pos file
% open file
fidPos = fopen([dataDir setNameNoExt '.pos'],'r','ieee-be'); %load pos
%'dirty' header read ( just for data start marker)
headerText = fread(fidPos,800,'int8'); %'dirty' header read
headerText = char(abs(headerText))';
ds = strfind(headerText,'data_start') + 10; % data start marker
% grab header in a more convenient format
frewind(fidPos);
headerPos = textscan(fidPos,'%s %[^\n]',27); % header has 28 lines - last is data start
headerPos = horzcat(headerPos{:});
% change n pos samples field
sRateInd = strcmp('sample_rate',headerPos(:,1)); % number of samples marker
sRatePos = sscanf(headerPos{sRateInd,2},'%d'); % get number of samples
% make index of valid samples
nPosSamplesInd = false(1,trialDurOrg * sRatePos);
nPosTime = 1:sRatePos * trialDurOrg;
for a = 1:size(timeS,1)
    nPosSamplesInd(nPosTime > timeS(a,1)*50 & nPosTime <= timeS(a,2)*50) = true;
end
nPosSamplesNew = sum(nPosSamplesInd);

% read pos data
fseek(fidPos,ds-1,'bof'); % Set file position indicator at data start
tempPosData = fread(fidPos,'int16');
tempPosData = tempPosData(1:trialDurOrg * sRatePos * 10);

% truncate pos data %
% rehape pos sample index
nPosSamplesInd = repmat(nPosSamplesInd,10,1); %each pos sample consists of 10 samples - 1 x 4 byte framecount + 8 x 2byte samples
nPosSamplesInd = nPosSamplesInd(:);
tempPosData = reshape( tempPosData(nPosSamplesInd),10,nPosSamplesNew); % format is data x n pos samples
tempPosData = tempPosData(3:10,:); % ignore frame count
frameCount = int32(1:nPosSamplesNew); % make new frame count
fclose(fidPos);

% back up original file
movefile([dataDir setNameNoExt '.pos'],[backUpDir setNameNoExt '.pos']); 

% change a few header fields
durInd = strcmp('duration',headerPos(:,1));
headerPos{durInd,2} = num2str(trialDurNew);
nSamplesInd = strcmp('num_pos_samples',headerPos(:,1)); % number of samples marker
headerPos{nSamplesInd,2} = num2str(nPosSamplesNew);
timeInd = strcmp(headerPos(:,1),'trial_time');
headerPos{timeInd,2} = char(trialTimeNew);
% write pos file - header first
fid = fopen([dataDir setNameNoExt '.pos'],'w');

for head = 1:length(headerPos)
    fprintf(fid,'%s %s\r\n',headerPos{head,1},headerPos{head,2});
end
fprintf(fid,'%s','data_start'); %add data start marker
% get data in the right format 
posData = int16(tempPosData); 
[~,~,endian] = computer; % prob have to swap byte order here as win uses le by default
if strcmp(endian,'L')
    frameCount = swapbytes(frameCount);
    posData = swapbytes(posData);   
end
posDataWrite = reshape( typecast(frameCount,'int8'), 4, []); % cast as int8
posDataWrite(5:20,:) = reshape( typecast(posData(:),'int8'), 16, []); % cast as int8
% write
fwrite(fid,posDataWrite(:),'int8',0,'ieee-be'); %write data

fprintf(fid,'\r\n%s\r\n','data_end'); %add data end marker
fclose(fid);
fprintf('Wrote .pos file for truncated data: %s%s\n##########\n',dataDir, [setNameNoExt '.pos']);
clear timeInd ds nSamplesInd nPosSamplesNew fidPos fid head headerText sRateInd sRatePos durInd posDataWrite nPosTime 
          
   
%% truncate tet data
% grab all tetfiles
temp = dir([dataDir setNameNoExt '.*']);
allFiles = {temp(:).name};
tetFileInd = ~cellfun('isempty',regexp(allFiles,'(?<=[.])\d'));
tetfiles = allFiles(tetFileInd);
clear temp

indTime = cell(length(tetfiles),1);
for l = 1:length(tetfiles)
    %CB code - load
    [headerTet, timestampOrg, waveforms] = read_tetrode_file([dataDir tetfiles{l}]);
    
    indTime{l} = false(length(timestampOrg),1);
    timestamp = [];
    timeStampOffset = 0; % this offset will make sure that relative positions of time stamps remain unchanged in trial
    for a = 1:size(timeS,1)
        
        currTimeWindowInd = timestampOrg(:,1) >= timeS(a,1) & timestampOrg(:,1) <= timeS(a,2); % index for spike time stamps in current time window
        % need to create new time stamp to account for new length of
        % truncated data. We also need to make sure that we keep the
        % relative timings coherent (e.g. relationship to EEG)
        if sum(currTimeWindowInd) == 0
            % no spikes in current window
            timeStampOffset = timeStampOffset + timeS(a,2) - timeS(a,1);
        else
            % 3 cases when spikes are in current window
            if a == 1 && timeS(a,1) == 0
                % 1: first window + window starts with trial start
                timestamp = timestampOrg(currTimeWindowInd,1); % start == trial start
                timeStampOffset = timeS(a,2);
            elseif a == 1 && timeS(a,1) ~= 0
                % 2: first window + window starts not with trial start
                timestamp = timestampOrg(currTimeWindowInd,1) - timeS(a,1);  
                timeStampOffset = timeS(a,2) - timeS(a,1);
            else
                % 3: all other cases
                timestamp = [timestamp;  timeStampOffset + timestampOrg(currTimeWindowInd) - timeS(a,1)];
                timeStampOffset = timeStampOffset + timeS(a,2) - timeS(a,1); 
            end
            
        end
        % make index to filter waveforms (and keep it in case cut files
        % have to be truncated too); update every iteration
        indTime{l} = indTime{l} | currTimeWindowInd;
    end
    % restore convenient time stamp format
    timestamp = repmat(timestamp,1,4); 
    % truncate waveforms
    waveforms = waveforms(indTime{l},:,:);
    
    %make header - edit 3 fields
    durInd = strcmp(headerTet(:,1),'duration');
    headerTet{durInd,2} = num2str(trialDurNew);
    spkInd = strcmp(headerTet(:,1),'num_spikes');
    headerTet{spkInd,2} = num2str(size(waveforms,1));
    timeInd = strcmp(headerTet(:,1),'trial_time');
    headerTet{timeInd,2} = char(trialTimeNew);
    
    %write truncated data
    %grab time base
    tBaseInd = strcmp(headerTet(:,1),'timebase');
    tBase = sscanf(headerTet{tBaseInd,2},'%d');
    %get wf's into right shape
    waveforms = reshape(waveforms,size(waveforms,1),200);
    %retransform time stamp into 4 byte
    timestamp = int32(timestamp(:,1) * tBase);
    %needs to be big endian format
    if strcmp(endian,'L')
        timestamp = swapbytes(timestamp);
    end
    timestamp = reshape( typecast(timestamp,'int8'),4,[]);
    %assign into output
    tempData = int8(zeros(216,size(waveforms,1)));
    tempData([5:54, 59:108, 113:162, 167:216],:) = waveforms';
    tempData([1:4, 55:58, 109:112, 163:166],:) = repmat(timestamp,[4,1]);
    
%     % back up original file
    movefile([dataDir tetfiles{l}],[backUpDir tetfiles{l}]); 
    
    %first write new header
    fid = fopen([dataDir tetfiles{l}],'w');
    for i = 1:length(headerTet)
        fprintf(fid,'%s %s\r\n',headerTet{i,1},headerTet{i,2});
    end
    fprintf(fid,'%s','data_start'); %add data marker
    fwrite(fid,tempData(:),'int8',0,'ieee-be'); %write data
    fprintf(fid,'\r\n%s\r\n','data_end'); %add data end marker
    fclose(fid);
    
    fprintf('Wrote tet file for truncated data: %s%s\n##########\n',dataDir, tetfiles{l});

end 
clear tetFileInd fid tBase endian spkInd timeInd durInd tetFileInd temp allFiles timestampOrg

%% truncate eeg/egf files
% grab all eeg/egf iles
temp = dir([dataDir setNameNoExt '.eeg*']);
eegFiles = {temp(:).name};
temp = dir([dataDir setNameNoExt '.egf*']);
egfFiles = {temp(:).name};

for n = 1:size(eegFiles,2)   
    % do EEGs
    fidEEG = fopen([dataDir eegFiles{n}],'r','ieee-be'); %open EEG file
    binData = fread(fidEEG,'int8'); % binary data
    ds = strfind(binData','data_start') + 10; % data start marker
    de = strfind(binData','data_end') - 3; % data end marker
    eegData = binData(ds:de); % grab voltage data

    % grab header and edit a few lines- only need to do once
    if n == 1
        frewind(fidEEG);
        headerEEG = textscan(fidEEG,'%s %[^\n]',11); % header has 12 lines - last is data start
        headerEEG = horzcat(headerEEG{:});
        sRateInd = strcmp('sample_rate',headerEEG(:,1)); % number of samples marker
        sRateEEG = sscanf(headerEEG{sRateInd,2},'%d'); % get number of samples
        % make index of valid samples
        nEEGSamplesInd = false(1,trialDurOrg * sRateEEG);
        nEEGTime = 1:sRateEEG*trialDurOrg;
        for a = 1:size(timeS,1)
            nEEGSamplesInd(nEEGTime > timeS(a,1) * sRateEEG & nEEGTime <= timeS(a,2) * sRateEEG) = true;
        end
        nEEGSamplesNew = sum(nEEGSamplesInd);
        nSamplesInd = strcmp('num_EEG_samples',headerEEG(:,1)); % number of samples marker
        headerEEG{nSamplesInd,2}= num2str(nEEGSamplesNew);       
        durInd = strcmp('duration',headerEEG(:,1));
        headerEEG{durInd,2} = num2str(trialDurNew);
        timeInd = strcmp(headerEEG(:,1),'trial_time');
        headerEEG{timeInd,2} = char(trialTimeNew);
    end
    eegData = eegData(nEEGSamplesInd); % truncate data

    fclose(fidEEG);
    
    % back up original file
    movefile([dataDir eegFiles{n}],[backUpDir eegFiles{n}]); 
    
    %do egf files if there are any - pretty much same thing
    if ~isempty(egfFiles)
        fidEGF = fopen([dataDir egfFiles{n}],'r','ieee-be'); %open header
        % get data start marker
        binDataEGF = fread(fidEGF,'int8');
        ds = strfind(binDataEGF','data_start') + 10; % data start marker
        %do header
        if n == 1
            frewind(fidEGF);
            headerEGF = textscan(fidEGF,'%s %[^\n]',10);
            headerEGF = horzcat(headerEGF{:});
            sRateInd = strcmp('sample_rate',headerEGF(:,1)); % number of samples marker
            sRateEGF = sscanf(headerEGF{sRateInd,2},'%d'); % get number of samples
            % make index of valid samples
            nEGFSamplesInd = false(1,trialDurOrg * sRateEGF);
            nEGFTime = 1:sRateEGF*trialDurOrg;
            for a = 1:size(timeS,1)
                nEGFSamplesInd(nEGFTime > timeS(a,1) * sRateEGF & nEGFTime <= timeS(a,2) * sRateEGF) = true;
            end
            nEGFSamplesNew = sum(nEGFSamplesInd);
            nSamplesInd = strcmp('num_EGF_samples',headerEGF(:,1)); % number of samples marker
            headerEGF{nSamplesInd,2} = num2str(nEGFSamplesNew);
            durInd = strcmp('duration',headerEGF(:,1));
            headerEGF{durInd,2} = num2str(trialDurNew);
            timeInd = strcmp(headerEGF(:,1),'trial_time');
            headerEGF{timeInd,2} = char(trialTimeNew);
        end
        
        % grab actual voltage data
        fseek(fidEGF,ds-1,'bof');
        egfData = fread(fidEGF,'int16'); %re-read as int16
        egfData = egfData(nEGFSamplesInd); % truncate data

        fclose(fidEGF);
        
        % back up original file
        movefile([dataDir egfFiles{n}],[backUpDir egfFiles{n}]); 
    end
    
    % write eeg file
    fid = fopen([dataDir eegFiles{n}],'w');
    for head = 1:length(headerEEG)
        fprintf(fid,'%s %s\r\n',headerEEG{head,1},headerEEG{head,2});
    end
    fprintf(fid,'%s','data_start'); %add data start marker
    fwrite(fid,eegData,'int8',0,'ieee-be'); %write data
    fprintf(fid,'\r\n%s\r\n','data_end'); %add data end marker
    fclose(fid);
    fprintf('Wrote .eeg file for truncated data: %s%s\n##########\n',dataDir, eegFiles{n});
  
    % write egf (if there are any)
    if ~isempty(egfFiles)
        fid = fopen([dataDir egfFiles{n}],'w');
        for head = 1:length(headerEGF)
            fprintf(fid,'%s %s\r\n',headerEGF{head,1},headerEGF{head,2});
        end
        fprintf(fid,'%s','data_start'); %add data start marker
        fwrite(fid,egfData,'int16',0,'ieee-be'); %write data
        fprintf(fid,'\r\n%s\r\n','data_end'); %add data end marker
        fclose(fid);
        fprintf('Wrote .egf file for truncated data: %s%s\n##########\n',dataDir, egfFiles{n});
    end
end

clear durInd eegNum nEGFSamplesNew nEEgSamplesNew nSamplesInd eegFiles egfFiles ds de timeInd fidEEG fidEGF fid nSamplesInd head sRateInd nEEGTime ...
      nEGFTime

%% truncate cut files
% easiest to use UI due to possibility of multiple cut files for
% datasets;
if cutFileFlag
    [cutNames, dummy] = uigetfile([dataDir setNameNoExt '*.cut'],'Select .cut Files of Trials to Join','MultiSelect','on');
    if dummy == 0
        warning('No cut files selected. Your journey ends here my friend.');
        fclose all;
        return
    end
    
    % some sanity check
    if length(cutNames) ~= length(tetfiles)
        warning('You selected a different number of cut files and tetrode files! Truncating cutfiles aborted');
        return
    end
    
    for p = 1:length(cutNames)
        fid = fopen([dataDir cutNames{p}]); %open
        headerCut = textscan(fid,'%s','delimiter','\n'); %read string
        %last line of header
        headerEnd = find(strncmp(headerCut{1},'Exact_cut_for',13));
        %header string - keep for later
        headerStr = headerCut{1}(1:headerEnd-1);
        %read cluster IDs, skip header
        frewind(fid);
        cluIDs = textscan(fid,'%n','headerlines',headerEnd,'delimiter','\b');
        cluIDs = cluIDs(indTime{p}); %truncated
        
        fclose(fid);
        
        % back up original file
        movefile([dataDir cutNames{p}],[backUpDir cutNames{p}]);
        
        headerStr = [headerStr(:);['Exact_cut_for: ' dataDir cutNames{p} ' spikes: ' num2str(size(waveforms,1))]];
        cutFile = [headerStr(:); cellstr(num2str(cluIDs))];
        % write file
        fid = fopen([dataDir cutNames{p}],'w'); %open file for writing
        headerLength = length(headerStr);
        
        % write file
        for r = 1:length(cutFile)
            %header lines with 'min' or 'max' values for centres - need to add
            %white space
            if ~isempty(strfind(cutFile{r},'min')) || ~isempty(strfind(cutFile{r},'max'))
                fprintf(fid,'              %s\r\n',cutFile{r});
                %cluster IDs
            elseif r > length(headerStr)
                % write cluster data here
                if ~mod(r-headerLength,25) %there is a carriage return every 25 entries
                    fprintf(fid,'%s\r\n',cutFile{r});
                else
                    fprintf(fid,'%s ',cutFile{r});
                end
                %any other header line
            else
                fprintf(fid,'%s\r\n',cutFile{r});
            end
        end
        fclose(fid);
        fprintf('Wrote .cut file for truncated data: %s%s\n##########\n',dataDir, cutNames{p});
    end    
end
fclose all;
end



