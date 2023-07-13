function splitCutFileFromMaster( mode,tetrodes,varargin )
% Splitting cluster files - improved version of Shaz's original function.
% Will make individual (i.e. for each trial) cut files from a master cut
% We assume that trial data names are of format:
% dateStamp(american format) - trial index (letter) - underscore - .....
% Note: this is slow when reading/writing from server
%
% Usage:    splitCutFileFromMaster( mode )
%           splitCutFileFromMaster(mode, tetrodes, optionalInputStruct)
%           splitCutFileFromMaster(mode, [], optionalInputStruct)
%           splitCutFileFromMaster(mode, tetrodes, 'inputName', inputVal, .. etc .. )
%           splitCutFileFromMaster(mode, [], 'inputName', inputVal, .. etc .. )
%
% Inputs:   mode - 'hc' = HardCoded (use info in params) or 'ui'= UserInput
%                  (open user dialog to select files)
%           tetrodes - numeric array of tetrodes to split (any in [1 2 3 4
%                      5 6 7 8 9 10 11 12 13 14 15 16]); if left empty/not supplied then we will
%                      assume the full set of available tetrodes 
% - THE WAY YOU RUN THIS IS:   splitCutFileFromMaster('hc',[1 2 3 4 5 6 7 8 9 10 11 12]) IV  20/06/23
%
% Optional inputs/analysis parameters (supply as " ,'fieldname',value, " comma-separated list) :
%
%           prms.dir = path to master directory
%           prms.Tnames = cell array with 1st & last trial in sequence
%           prms.TNames2exclude = cell array with trial names to be excluded   
%           prms.tets = [1 2 3 4 5 6 7 8]; % n of tetrodes; default set
%           prms.masterCutID = string to identify master cut ID
%           prms.cutFileExt = string to be added to cut file @ .._(tetNum)xx.cut 
%           prms.addStr = string to be added to cut file @ ..xx_(tetNum)...
%
% LM 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
% specify directory with .set files
prms.dir = 'S:\DBIO_TKFC_SKGTIVA\thesis_data\tetrode_data\r804\';%' %Z:\lmuessig\!postDoc\recording_data\CA1\; 
% specify names of trial files
prms.Tnames = {'190528a_famBox','190528d_famBox'}; % first and last trial in sequence _famBox _sleepHP _sqTrack 
prms.TNames2exclude = {''};%{'201218f_sleepHP'}; % trials that should be excluded
% specify tetrodes  
prms.tets = [5 6 7 8 9 10 11 12]; % if empty, will do all files, otherwise can also specify numbers as array (needs to be row vector)
% string which identifies master cut files
prms.masterCutID = '_alltrials_'; %'allTrials_sqTrack'; 'allLTM; 'allReact'; 'allTrials_novLinTrack'; 'allLinear';'DGCA3';
prms.cutAddStrMaster = []; %'fin'; 
prms.cutFileExt =  '';%'firstcut';%extra str after tet Number in cutfile - leave blank if your naming system is different
% added to cut file name (after ..._(tetN)
prms.addStr =  '_DGCA3';%''; '_sqTrack'; '_LTM'; '_react';'_novLinTrack'; '_phaseShift';

% parse optional inputs
if nargin>2 && isstruct(varargin{1})
    optIn = varargin{1};
    f=fieldnames(optIn);
    for i=1:length(f);   prms.(f{i}) = optIn.(f{i});   end
elseif nargin>2 && ischar(varargin{1})
    for i=1:2:length(varargin)
        prms.(varargin{i}) = varargin{i+1};
    end
end

%% parsing filenames 
if strcmp(mode,'hc')
    %grab date stamp from trial name - assuming YYMMDD
    dateStamp = prms.Tnames{1}(1:6);
    %retrieve filenames to load
    files = dir([prms.dir '*' dateStamp '*.set']); %only check for relevant files, i.e .set files containing the date stamp
    startInd = ~cellfun('isempty',regexp({files.name},['\w*' prms.Tnames{1} '\w*'])); %find first file of trial sequence
    startInd = find(startInd); %numeric index
    endInd = ~cellfun('isempty',regexp({files.name},['\w*' prms.Tnames{2} '\w*'])); %find last file of trial sequence
    endInd = find(endInd);
    files = files(startInd:endInd);
    files = {files(:).name};
elseif strcmp(mode,'ui')
    [files, prms.dir] = uigetfile('*.set','Select .set Files of Trials','MultiSelect','on');
    if ~ischar(prms.dir)
        error('Something went wrong with the file selection! Can''t do anything man.');
    end
    prms.TNames2exclude = {}; % should not be used when loading from UI 
    dateStamp = files{1}(1:6); %grab date stamp from trial name - assuming YYMMDD
end
% remove extension
for j = 1:length(files)
    [~,Tnames{j},~] = fileparts(files{j}); %remove extension
end

%remove trials that should be excluded
ind = ~ismember(Tnames,prms.TNames2exclude);
Tnames = Tnames(ind);

if nargin == 1 || isempty(tetrodes)
    [~,tetrodes] = findDACQFiles( prms.dir, Tnames{1}, 'tet' );
end

%pre-allocate
samplesTrial = zeros(length(Tnames),1);

%% loop through data (tetrodes)

for i = 1:length(tetrodes)
    
    %grab  
    cutf = dir([prms.dir '*' dateStamp '*_' num2str(tetrodes(i)) prms.cutAddStrMaster '.cut']);
    %find index for master cut file
    masterInd = ~cellfun('isempty',strfind({cutf(:).name},prms.masterCutID));

    %skip if no master cut for current tetrode
    if sum(masterInd) == 0
        fprintf('No master cut file available for tetrode %d\n',tetrodes(i));
        continue;        
    end
    %load master cut
    masterCutFile = [prms.dir cutf(masterInd).name];
    fid = fopen(masterCutFile); %open
    header = textscan(fid,'%s','delimiter','\n'); %read string
%     fclose(fid);
    %last line of header
    headerEnd = find(strncmp(header{1},'Exact_cut_for',13));
    %header string - keep for later
    headerStr = header{1}(1:headerEnd-1);
%     fid = fopen(masterCutFile); %open again
    frewind(fid);
    %read cluster IDs, skip header
    clall = textscan(fid,'%n','headerlines',headerEnd,'delimiter','\b');
    clall = [clall{:}];
    fclose(fid);

    %grab number of time stamps from tetrode files
    for t = 1:length(Tnames)
        wavs = load_spikes([prms.dir Tnames{t} '.' num2str(tetrodes(i))]);
        samplesTrial(t) = length(wavs.times);
    end

    %make cut file indices
    cutind = [0;cumsum(samplesTrial)];
    %write individual cluster files derived from file above...
    for t = 1:length(Tnames)
        
        i1 = cutind(t) + 1;
        i2 = cutind(t+1);
        %get cluster IDs for current trial
        clTrial = clall(i1:i2);
        %output name
        TrialcutfileName = [prms.dir Tnames{t} prms.addStr '_' num2str(tetrodes(i)) prms.cutFileExt '.cut'];
        %add last line for header in current trial
        headerStrTrial = [headerStr(:);['Exact_cut_for: ' Tnames{t} ' spikes: ' num2str(samplesTrial(t))]];
        %write file
        fName = writeCutFile( clTrial,headerStrTrial,TrialcutfileName );
        
        fprintf('Radikalisnky! Cut file successfully generated: %s\n',fName);
    end
end
end



