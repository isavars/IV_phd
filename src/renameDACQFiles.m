function renameDACQFiles
% rename all DACQ files that belong to a specific trial (in case of typo/error
% when saving data
%
% LM 2018 


defaultDir = 'Z:\IsabellaV\recording_data\'; %'Z:\lmuessig\!postDoc\recording_data\CA1\'; 'Z:\lmuessig\!postDoc\recording_data\CA1\'; 'D:\tempRecordingData\'; 'D:\tempAdultData\';

[fName, dataDir] = uigetfile([defaultDir '*.*'],'Select File From Group That You Want Renamed');
if ~ischar(dataDir)
    warning('Loading was cancelled. Can''t continue. Nooooo.');
    return
end
% get all files from group
[~,fNameNoExt,ext] = fileparts(fName);
if strcmp(ext,'.cut')
    warning('Don''t use cut files to get file name mate. Any other file types are fine, yeah!');
    return
end
fList = dir([dataDir fNameNoExt '.*']);
% new name (Be careful not to have two files with the same name while remaning. They will be overwritten with no warning)
newFName = inputdlg('New file name: ','Please fill in',[1 50],{fNameNoExt});

% for backup of original data
backupDir = [dataDir 'tempBackup\'];
if ~isdir(backupDir)
    mkdir(backupDir);
end
% rename and move
for i = 1:length(fList)
    [~,~,ext] = fileparts(fList(i).name);
    copyfile([dataDir fList(i).name],[dataDir newFName{:} ext]);
    movefile([dataDir fList(i).name],[backupDir fList(i).name]);
end


end
