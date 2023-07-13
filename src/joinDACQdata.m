  function joinDACQdata(varargin)
% This function will concatonate all DACQ data across multiple trials into 
% joint files. Obviously this really only makes sense if the trials were 
% recorded in quick succession. 
% We assume that trial names are of general format:
% dateStamp(american format) - trial index (letter) - underscore - .....
% The output file will then be called 'dateStamp - prms.addStr - .ext
% We will load all file types from raw, create joint files and write those 
% to disk. 
%
% Note1: When the cut files are joined additional UI is needed, so hang in
% there and don't go anywhere.
%
% Note2: When selecting the .set files of trials to join, make sure they
% are in right (i.e. alphabetical order) as naming of joint files will 
% depend on that.
%
%  File formats:
%  Set - just Ascii header
%  Pos - header - data_start - samples(4byte time stamp + 8x2byte values) -
%        CR-LF-data_end-CR-LF
%  tet - header - data_start - samples(4byte time stamp-ch1-4byte time 
%        stamp-ch2...ch4) - CR-LF-data_end-CR-LF
%  eeg - header - data_start - samples (1byte) - CR-LF-data_end-CR-LF
%  egf - header - data_start - samples (2byte) - CR-LF-data_end-CR-LF
%  cut - header - data_start - samples (text)
%
%  Usage:  joinDACQdata
%          joinDACQdata( optionalInputStruct )
%          joinDACQdata( 'inputName', inputVal, .. etc .. )
%
%  Optional inputs/analysis parameters (supply as " ,'fieldname',value, " comma-separated list) :
%
%          prms.nTets = 1:8; %1:8 - tetrodes do do  
%          % default dir where files will be assumed to be; after 1st UI
%          % call dir will be switched (if different from default)
%          prms.defaultDir = 'C:\Users\LM\Dropbox\temp\New folder\';
%          prms.addTrialType = 1; %y/n flag
%          prms.cutTag1 = '_sqTrack';%'_sqTrack'; '';
%          prms.cutTag2 = 'fin'; % 'fin'; '';
%          
%
%  LM 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TO DO: 
% 1)add joining of clu files 
% 2)if settings are inconsistent across trials but only concern e.g. 1 EEG
% channel, skip joining for that channel only 

%% params
prms.loading = 'ui'; % 'ui' - user input; 'hc' - hard coded in prms.files
prms.files = '';
prms.nTets = []; % if empty, will do all files, otherwise can also specify numbers as array (needs to be row vector)
prms.defaultDir = 'Z:\IsabellaV\recording_data\';%'Z:\lmuessig\!postDoc\recording_data\CA1\'; %'Z:\lmuessig\!postDoc\recording_data\CA1\'; 'D:\tempRecordingData\';
% for naming stuff (most is idiosyncratic to my data, so simply adjust to your needs):
% for trial names
prms.addStr = 1; %y/n flag
%for cut files
prms.cutFlag = 0;
prms.cutTag1 = '_DGCA3'; %'_sqTrack'; % '_sqTrack'; '_react'; ''; - this is probably only useful for me
prms.cutTag2 = ''; % 'fin'; ''; - this is probably only useful for me

prms.leaveOverhang = 0;

%parse optional inputs
if nargin>0 && isstruct(varargin{1})
    optIn = varargin{1};
    f=fieldnames(optIn);
    for i=1:length(f);   prms.(f{i}) = optIn.(f{i});   end
elseif nargin>0 && ischar(varargin{1})
    for i=1:2:length(varargin)
        prms.(varargin{i}) = varargin{i+1};
    end
end

%% parse trial names etc
%parse trial names
addStr = '';
switch prms.loading
    case 'ui'
        [setNames, dataDir] = uigetfile([prms.defaultDir '*.set'],'Select .set Files of Trials to Join','MultiSelect','on'); 
        if ~ischar(dataDir)
            warning('Loading was cancelled. Can''t continue. Nooooo.');
            return
        end
        fNames = cell(length(setNames),1);
        for a = 1:length(setNames)
            [~,fNames{a},~] = fileparts([dataDir setNames{a}]); %remove extension
            strInd = strfind(fNames{a},'_');
            addStr = [addStr fNames{a}(strInd(1)-1)];
        end
    case 'hc'
        fNames = cell(length(prms.files),1);
        % get input in same format as if selected from UI
        for a = 1:length(prms.files)
            [dataDir,fNames{a},] = fileparts(prms.files{a});
            strInd = strfind(fNames{a},'_');
            addStr = [addStr fNames{a}(strInd(1)-1)];
        end
        dataDir = [dataDir '\'];
end
% grab trialType
trialType = fNames{1}(strInd(1)+1:end);

% output dir
writeDir = [dataDir]; % 'joinedFiles\'];
if ~isdir(writeDir)
    mkdir(writeDir);
end
% grab datestamp
dateStamp = fNames{1}(1:6);
% this might be idiosyncratic to my data, so skip if nameing convention is
% different
if prms.addStr
    % string to add to trial names
    addStr = [addStr '_' trialType];
end

flagOverhang = false(1,length(fNames)); % flag for trials with DACQ overhang - we wanna remove that from all data

fprintf('Concatonating DACQ data for trials %s-%s. ######\nJoined data will be put in directory %s ######\nand named %s....\n##########\n',fNames{1},fNames{end},writeDir,[dateStamp addStr]);

clear strInd setNames a

%% make joint set file
trialDur = nan(length(fNames),1);
tempJoinedSet = {};
for i = 1:length(fNames)
    [~, sFileText] = load_set([dataDir fNames{i}]); %load set files (use Scan function for convenience)
    durInd = strcmp(sFileText(:,1),'duration');
    temp = sscanf(sFileText{durInd,2},'%d');
    % deal with overhang in case there is some
    if ~rem(temp,10) || prms.leaveOverhang
        trialDur(i) = temp;
    else
        trialDur(i) = temp - rem(temp,10); % subtract overhang
        flagOverhang(i) = true; % flag it for other files
    end
    % take  both columns in first trial
    if i==1
        tempJoinedSet = cat(2,tempJoinedSet,sFileText);
    else
        tempJoinedSet = cat(2,tempJoinedSet,sFileText(:,2));
    end
end
% change a few fields in header
tempJoinedSet(durInd,2:end) = cellstr(num2str(sum(trialDur)));
timeInd = strcmp(tempJoinedSet(:,1),'trial_time');
trialTime = tempJoinedSet(timeInd,2);
tempJoinedSet(timeInd,2:end) = trialTime;
lastTrialInd = strcmp(tempJoinedSet(:,1),'lasttrialdatetime');
tempJoinedSet(lastTrialInd,2:end) = tempJoinedSet(lastTrialInd,2);
% check for inconsistencies across set files - these would indicate that
% some settings changed across trials. Maybe should return error here?
[~,~,checkNums] = unique(tempJoinedSet(:,2:end),'stable');  % trick: use equal numbers instead of strings
checkNums = reshape(checkNums,length(tempJoinedSet),size(tempJoinedSet,2)-1); % orig. format
checkSet = ~all(bsxfun(@eq,checkNums,checkNums(:,1)),2); % check
% print inconsistent fields
if sum(checkSet) ~= 0
    warning('Some properties across trials are not consistent! Better go and check that!');
    inconsistInd = find(checkSet==1);
    for in = 1:length(inconsistInd)
        fprintf('Field ''%s'' is inconsistent across trials.\n##########\n',tempJoinedSet{inconsistInd(in),1})
    end
end
setFinal = tempJoinedSet(:,1:2); %final set file
% write file
% make filename first
setFileName = [dateStamp addStr '.set'];
fid = fopen([writeDir setFileName],'w');
for j = 1:length(setFinal)
    fprintf(fid,'%s %s\r\n',setFinal{j,1},setFinal{j,2});
end
fclose(fid);
fprintf('Finished concatonating .set file for %s\n##########\n',[dateStamp addStr]);
clear setFileName tempJoinedSet timeInd lastTrialInd checkSet setFinal checkNums sFileText temp fid

%% concat spike data 
% use separate function
if isempty(prms.nTets)
    [~,Tet_iterator] = findDACQFiles( dataDir, fNames{1}, 'tet' );
else
    Tet_iterator = prms.nTets;
end
for k = Tet_iterator
    concat_TetrodeTrialData( k ,'tint','dataRootDir',dataDir,'fNames',fNames,'writeDir',writeDir,'addStr',addStr,'leaveOverhang',prms.leaveOverhang);
end
fprintf('Finished concatonating tetrode data for %s\n##########\n',[dateStamp addStr]);

%% pos file
% posFiles = cell(length(fNames),1);

posData = [];
frameCount = [];
nPosSamples = nan(length(fNames),1);
frameCountAdd = 0; % offset for frame count for trials > 1 
for l = 1:length(fNames)
%     temp = dir([dataDir fNames{l} '.pos']);
%     posFiles{l} = temp.name;
    % open file
    posFile = findDACQFiles( dataDir, fNames{l}, 'pos' );
    fidPos = fopen([dataDir posFile],'r','ieee-be'); %load pos 
    %'dirty' header read ( just for data start marker)
    headerText = fread(fidPos,800,'int8'); %'dirty' header read
    headerText = char(abs(headerText))';
    ds = strfind(headerText,'data_start') + 10; % data start marker
    % grab header in a more convenient format
    frewind(fidPos);
    headerPos = textscan(fidPos,'%s %[^\n]',27); % header has 28 lines - last is data start
    headerPos = horzcat(headerPos{:});
    % grab n pos samples field
    nSamplesInd = strcmp('num_pos_samples',headerPos(:,1)); % number of samples marker
    nPosSamples(l) = sscanf(headerPos{nSamplesInd,2},'%d'); % get number of samples
    % read pos data
    fseek(fidPos,ds-1,'bof'); % Set file position indicator at data start 
    tempPosData = fread(fidPos,'int16'); 
    % read data: 1 x 4 byte framecount +  8 x 2byte samples
    tempPosData = reshape( tempPosData(1:nPosSamples(l)*10),10,nPosSamples(l)); % format is data x n pos samples  
    tempPosData = tempPosData(3:10,:); % ignore frame count
    % read 4 byte frame counter - need to be treated separately
    fseek(fidPos,ds-1,'bof'); % go back 
    % read frame count and add appropriate amount of samples
    tempframeCount = fread(fidPos,nPosSamples(l),'int32',16) + frameCountAdd; 
    % remove DACQ overhang
    if flagOverhang(l)
        sRateInd = strcmp('sample_rate',headerPos(:,1)); % sample rate
        sRate = sscanf(headerPos{sRateInd,2},'%d');
        tempframeCount = tempframeCount(1:trialDur(l)*sRate,1); % truncate data
        tempPosData = tempPosData(:,1:trialDur(l)*sRate); % truncate data
        nPosSamples(l) = trialDur(l)*sRate; % update
    end
    % keep separate
    posData = [posData, tempPosData]; 
    frameCount = [frameCount; tempframeCount];    
    frameCountAdd = frameCount(end) + 1; 

    fclose(fidPos);
end
% change a few header fields
durInd = strcmp('duration',headerPos(:,1));
headerPos{durInd,2} = num2str(sum(trialDur));
headerPos{nSamplesInd,2} = num2str(sum(nPosSamples));
timeInd = strcmp(headerPos(:,1),'trial_time');
headerPos(timeInd,2) = trialTime;
% write pos file - header first
posFileName = [dateStamp addStr '.pos'];
fid = fopen([writeDir posFileName],'w');

for head = 1:length(headerPos)
    fprintf(fid,'%s %s\r\n',headerPos{head,1},headerPos{head,2});
end
fprintf(fid,'%s','data_start'); %add data start marker
% get data in the right format 
frameCount = int32(frameCount);
posData = int16(posData); 

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
fprintf('Finished concatonating pos data for %s\n##########\n',[dateStamp addStr]);
clear timeInd ds nSamplesInd nPosSamples posFiles temp headerPos fidPos posData tempData posFileName fid head headerText tempPosData tempframeCount...
        frameCountAdd posDataWrite frameCount

%% eeg/egf
% this is a bit more involved due to number of potential files
% collect eeg/egf file names
eegFiles = {}; egfFiles = {};
for m = 1:length(fNames)
%     temp = dir([dataDir fNames{m} '.eeg*']);
%     eegFiles(m,:) = {temp(:).name};
%     temp = dir([dataDir fNames{m} '.egf*']);
%     egfFiles(m,:) = {temp(:).name};
      eegFiles(m,:) = findDACQFiles( dataDir, fNames{m}, 'eeg' );
      egfFiles(m,:) = findDACQFiles( dataDir, fNames{m}, 'egf' );
end

for n = 1:size(eegFiles,2)
    %pre-allocate
    eegData = []; egfData = [];
    nEEGSamples = nan(size(eegFiles,1),1);
    nSamplesEGF = nan(size(egfFiles,1),1);
    % eeg numerator (a bit ugly)
    strInd = strfind(eegFiles{1,n},'eeg');
    eegNum = eegFiles{1,n}(strInd+3:end); 
    
    for o = 1:size(eegFiles,1) 
        % do EEGs
        fidEEG = fopen([dataDir eegFiles{o,n}],'r','ieee-be'); %open EEG file
        % get data start marker
        binData = fread(fidEEG,'int8');
        ds = strfind(binData','data_start') + 10; % data start marker
        de = strfind(binData','data_end') - 3; % data start marker
        % grab header in a more convenient format
        frewind(fidEEG);
        headerEEG = textscan(fidEEG,'%s %[^\n]',11); % header has 12 lines - last is data start
        headerEEG = horzcat(headerEEG{:});
        nSamplesInd = strcmp('num_EEG_samples',headerEEG(:,1)); % number of samples marker
  
        %grab actual voltage data
        if flagOverhang(o)
            sRateInd = strcmp('sample_rate',headerEEG(:,1)); % number of samples marker
            sRate = sscanf(headerEEG{sRateInd,2},'%d');
            binData = binData(ds:de); % grab binary data
            eegData = [eegData; binData(1:trialDur(o)*sRate)]; % truncate data
            nEEGSamples(o) = trialDur(o)*sRate;
        else
            nEEGSamples(o) = sscanf(headerEEG{nSamplesInd,2},'%d'); % get number of samples
            eegData = [eegData; binData(ds:de)]; % grab binary data
        end
        fclose(fidEEG);
        
        %do egf files if there are any - pretty much same thing
        if ~isempty(egfFiles)
            fidEGF = fopen([dataDir egfFiles{o,n}],'r','ieee-be'); %open header
            % need to re-do as EGF header has one line less
            headerText = fread(fidEGF,'int8');
            ds = strfind(headerText','data_start') + 10; % data start marker 
            % more convenient format
            frewind(fidEGF);
            tempHeader = textscan(fidEGF,'%s %[^\n]',10);
            tempHeader = horzcat(tempHeader{:});
            nSamplesIndEGF = strcmp('num_EGF_samples',tempHeader(:,1)); 

            % grab actual voltage data
            fseek(fidEGF,ds-1,'bof');  
            tempDataEGF = fread(fidEGF,'int16'); %re-read as int16
            % deal with overhang
            if flagOverhang(o)
                sRateInd = strcmp('sample_rate',tempHeader(:,1)); % sample rate
                sRate = sscanf(tempHeader{sRateInd,2},'%d');
                nSamplesEGF(o) = trialDur(o)*sRate;
                tempDataEGF = tempDataEGF(1:nSamplesEGF(o));
                egfData = [egfData; tempDataEGF(1:trialDur(o)*sRate)];
            else
                nSamplesEGF(o) = sscanf(tempHeader{nSamplesIndEGF,2},'%d');
                egfData = [egfData; tempDataEGF(1:nSamplesEGF(o))];
            end
            fclose(fidEGF);
        end
        
    end
    % change a few header fields for EEG
    durInd = strcmp('duration',headerEEG(:,1));
    headerEEG{durInd,2} = num2str(sum(trialDur));
    headerEEG{nSamplesInd,2} = num2str(sum(nEEGSamples));
    timeInd = strcmp(headerEEG(:,1),'trial_time');
    headerEEG(timeInd,2) = trialTime;
    % write eeg file
    currEEGFileName = [dateStamp addStr '.eeg' eegNum];
    fid = fopen([writeDir currEEGFileName],'w');
    for head = 1:length(headerEEG)
        fprintf(fid,'%s %s\r\n',headerEEG{head,1},headerEEG{head,2});
    end
    fprintf(fid,'%s','data_start'); %add data start marker
    fwrite(fid,eegData,'int8',0,'ieee-be'); %write data
    fprintf(fid,'\r\n%s\r\n','data_end'); %add data end marker
    fclose(fid);
    % write egf (if there are any)
    if ~isempty(egfFiles)
        % write egf
        headerEGF = headerEEG; %headers are quite similar
        % remove 'EEG_samples_per_position' field as not part of egf header
        eegSamplesInd = strcmp('EEG_samples_per_position',headerEGF(:,1));
        headerEGF = headerEGF(~eegSamplesInd,:);
        bytesInd = strcmp('bytes_per_sample',headerEGF(:,1));
        headerEGF{bytesInd,2} = '2'; % just hard code
        sRateInd = strcmp('sample_rate',headerEGF(:,1));
        headerEGF{sRateInd,2} = '4800 Hz'; % just hard code
        headerEGF{nSamplesIndEGF,1} = 'num_EGF_samples';
        headerEGF{nSamplesIndEGF,2} = num2str(sum(nSamplesEGF));
        % write file
        currEGFFileName = [dateStamp addStr '.egf' eegNum];
        fid = fopen([writeDir currEGFFileName],'w');
        for head = 1:length(headerEGF)
            fprintf(fid,'%s %s\r\n',headerEGF{head,1},headerEGF{head,2});
        end
        fprintf(fid,'%s','data_start'); %add data start marker
        fwrite(fid,egfData,'int16',0,'ieee-be'); %write data
        fprintf(fid,'\r\n%s\r\n','data_end'); %add data end marker
        fclose(fid);
    end
end
fprintf('Finished concatonating eeg\\egf data for %s\n##########\n',[dateStamp addStr]);

clear headerEGF bytesInd headerEEG durInd currEEGFileName currEGFFileName eegNum egfData eegData nSamples nSamplesEGF tempDataEGF tempData nSamplesInd eegFiles egfFiles nSamplesIndEGF ds ...
    timeInd fidEEG fidEGF fid temp eegSamplesInd head sRateInd headerText nEEGSamples


%% join cut files
if prms.cutFlag
    % easysiest to use UI due to possibility of multiple cut files for
    % datasets; also do 1 trial at a time so one doesn't need to go through a
    % long list of cut files
    cutNames = {};
    for nFiles = 1:length(fNames)
        [temp, ~] = uigetfile([dataDir fNames{nFiles} '*.cut'],'Select .cut Files of Trials to Join','MultiSelect','on');
        cutNames = [cutNames,temp];
    end
    clear temp nFiles
    % some sanity check
    if length(cutNames) ~= Tet_iterator * length(fNames)
        error('You didn''t select the correct amount of cut files man. Damn!')
    end
    
    for p = Tet_iterator
        ind = ~cellfun('isempty',strfind(cutNames,['_' num2str(p)]));%%$%%%% NEEDS FIX FOR 16 TET DATA
        tempCutNames = cutNames(ind); % cut files for current tetrode
        
        cluIDs = [];
        nSpkTrial = nan(length(tempCutNames),1);
        for q = 1:length(fNames)
            fid = fopen([dataDir tempCutNames{q}]); %open
            headerCut = textscan(fid,'%s','delimiter','\n'); %read string
            %last line of header
            headerEnd = find(strncmp(headerCut{1},'Exact_cut_for',13));
            %header string - keep for later
            headerStr = headerCut{1}(1:headerEnd-1);
            %read cluster IDs, skip header
            frewind(fid);
            temp = textscan(fid,'%n','headerlines',headerEnd,'delimiter','\b');
            fclose(fid);
            
            % grab number of time stamps from tetrode files
            wavs = load_spikes([dataDir fNames{q} '.' num2str(p)]);
            % deal with the possibility of DACQ overhang
            if flagOverhang(q)
                ind = wavs.times <= trialDur(q);
                nSpkTrial(q) = length(wavs.times(ind));
                cluIDs = [cluIDs; temp{:}(ind)];
            else
                cluIDs = [cluIDs; temp{:}];
                nSpkTrial(q) = length(wavs.times);
            end
            
        end
        %final file
        currCutFileName = [dateStamp addStr prms.cutTag1 '_' num2str(p) prms.cutTag2 '.cut'];
        headerStrFin = [headerStr(:);['Exact_cut_for: ' currCutFileName ' spikes: ' num2str(sum(nSpkTrial))]];
        cutFile = [headerStrFin(:); cellstr(num2str(cluIDs))];
        % write file
        fid = fopen([writeDir currCutFileName],'w'); %open file for writing
        headerLength = length(headerStrFin);
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
    end
    fprintf('Finished concatonating .cut files for %s.\n##########\nALL done now!\n##########\n',[dateStamp addStr]);
end
fclose all;
end
