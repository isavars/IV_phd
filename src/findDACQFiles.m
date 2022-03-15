function [ fNames, validTrodes ] = findDACQFiles( dataDir, identifier, type, removeEXT )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% default - keep extension
if nargin == 3 || strcmp(type,'tet')
    removeEXT = 0;
end

validTrodes = [];

% find specific type
switch type
    case 'set'
        temp = dir([dataDir identifier '*.set']);
        fNames = {temp(:).name};   
    case 'pos'
        temp = dir([dataDir identifier '*.pos']);
        fNames = {temp(:).name};       
    case 'eeg'
        temp = dir([dataDir identifier '*.eeg*']);
        sz   = [temp.bytes];    % Exclude 'dummy' eeg files, which are just a partial text header.
        temp = temp( sz>1000 );  %
        fNames = {temp(:).name};   
    case 'egf'
        temp = dir([dataDir identifier '*.egf*']);
        sz   = [temp.bytes];
        temp = temp( sz>1000 );
        fNames = {temp(:).name};           
    case 'tet'
        % a bit more involved to find all tetrode files
        temp = dir([dataDir identifier '*']);
        allFiles = {temp(:).name};
        tetFileInd = ~cellfun('isempty',regexp(allFiles,'(?<=[.])\d')); 
        % clu files also have integer after a '.', so we want to exclude them.
        % TW edit: Laurenz, you forgot about fet, klg, fmask - you also need to exclude these.
        cluFilesInd = false( size(allFiles) );
        cluFileTags = {'.clu.', '.fet.', '.klg.', '.fmask.'};
        for ii=1:length(cluFileTags)
            cluFilesInd = cluFilesInd | ~cellfun('isempty',strfind(allFiles, cluFileTags{ii})); 
        end
        fNames = allFiles(tetFileInd & ~cluFilesInd); % all tetrode file 

        % can also get numeric list of tet numbers 
        if nargout == 2
            splitNames = regexp(fNames, '[.]', 'split'); % split into name + extension
            splitNames = vertcat(splitNames{:});
            validTrodes = sort(str2double(splitNames(:,2)))'; % tetrode numbers that were recorded;
        end  
    otherwise        
        error(['"' type '"' 'is not a valid file type! Embarrassing dude'])
end

if length(fNames) == 1
    fNames = fNames{:};
end

% remove file extension, if desired
if removeEXT
    splitNames = regexp(fNames, '[.]', 'split'); % split into name + extension
    splitNames = vertcat(splitNames{:});
    fNames = splitNames(:,1)'; 
end
end
