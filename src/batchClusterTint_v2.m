%batchClusterTint_v2
% script to batch cluster data using Dan's fancy KK algorithm - it will
% first concatonate terodes across trials and then call Dan's fancy KK for
% all data. The script assumes trial data is of format:
% dateStamp(american format) - trial index (letter) - underscore - .....
% The output file will then be called 'dateStamp - addStr - .(tetNum)
% The clu files will be called 'dateStamp - addStr - .(tetNum).clu (as
% per standard KK convention
% Note: Dan's KK doesn't work if there are special characters in any folder 
% name that belongs to the path
%
% LM 2016 

%% params
defaultDir = 'S:\DBIO_TKFC_SKGTIVA\thesis_data\tetrode_data\r804\';
writeDir = 'S:\DBIO_TKFC_SKGTIVA\thesis_data\tetrode_data\r804\';%'D:\tempKK\'; % data will be writen to writeDir\folderNames
nTetrodes = []; %if empty, will do all files, otherwise can also specify numbers as array (needs to be row vector)
fancyKKpath = 'C:\Program Files\Axona\fancykk\';%'D:\fancykk\'; %path to fancy KK folder on disk
TempFolderPath = 'C:\Users\Isabella\AppData\Local\Temp\'; %path to temp folder (usually hidden) in windows (updating username should do the trick)

% select data
[fNames, dataDir] = uigetfile([defaultDir '*.set'],'Select .set Files of Trials to Cluster','MultiSelect','on');
if ~ischar(dataDir)
    warning('Loading was canceled. Can''t continue. Nooooo.');
    return
end
% for output dir
ratID = inputdlg('Animal ID please: ','Please fill in',[1 25],{'r'});
tempWriteDir = [writeDir ratID{:} '\'];
if ~isdir(tempWriteDir)
    mkdir(tempWriteDir);
end
% for naming of master cut
addStr = inputdlg('Please enter string for master cut files: ','Master Cut File Identifier',[1 50]);

% you might need to change next line if your data is named following a
% different convention
dateStamp = fNames{1}(1:6); % assuming beginning of filename is YYMMDD
 
% remove extension
finFNames = cell(1,length(fNames));
for i = 1:length(fNames)
    [~,finFNames{i},~] = fileparts(fNames{i}); %remove extension
end
% some feedback to check that there is no mistake in trial selection
fprintf('These trials will be clustered boss:\n');
for j = 1:length(finFNames)
    fprintf('## %s ##\n',finFNames{j});
end

cd(tempWriteDir);

% gather tetrodes
if isempty(nTetrodes)
    [~,nTetrodes] = findDACQFiles( dataDir, finFNames{1}, 'tet' );
end

% concat tetrode data
fNameOut = cell(1,length(nTetrodes));
for k = nTetrodes
    fNameOut{k} = [ dateStamp addStr{:} '.' num2str(k) ];
    if exist([tempWriteDir fNameOut{k} ],'file') == 2
        
        if ~isdir([tempWriteDir '\old'])
            mkdir([tempWriteDir '\old']);
        end
        movefile(fNameOut{k},[tempWriteDir '\old\' fNameOut{k}]);
        warning([fNameOut{k} ' is exisiting file. Old file moved to ' tempWriteDir '\old\' fNameOut{k}]);
    end
    
    concat_TetrodeTrialData( k, 'tint','dataRootDir',dataDir,'fNames',finFNames,'writeDir',tempWriteDir,'addStr',addStr{:});
    
end

% run Dan's KK
if length(fNameOut) > 1
    % write list file for fancy KK call - needs to be stored in windows
    % %TEMP% folder
    fid = fopen([TempFolderPath 'fancy_kk_commands.list'],'w');
    for l = nTetrodes
        fprintf(fid,'"%s%s"\r\n',tempWriteDir,fNameOut{l});
    end
    fclose(fid);
    system([fancyKKpath 'runkk.bat'],'-echo'); %call fancy KK
else
    % in case there is only one tetrode to run, call directly
    system([fancyKKpath 'runkk.bat ' '"' fNameOut{:} '"'],'-echo'); %call fancy KK
end
% clean up workspace
clear tempWriteDir files finFNames nTetrodes fancyKKpath dateStamp i j k nTetrodes ratID dataDir fNames defaultDir 