function cutFileName = writeCutFile( clustIDs,headerStr,cutFileName )
% write a cut file from scratch by supplying the cluster ID and header. The
% file will have exactly same format as if written in TINT
%
% Usage: [ fName ] = writeCutFile( clustIDs,headerStr,fName )
%
% Inputs:   clustIDs    - n-by-1 array of cluster IDs
%           header      - cell array with header strings
%           cutFileName - filename of cut file - if path is not specified
%                         cd path will be added 
% Outputs:  cutFileName - path + filename of cut file that was generated
%
% LM 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make sure filename includes path to directory
% add path (to cd) to filename if path is not already included in 'cutFileName'
if isempty(strfind(cutFileName,'\'))
    cutFileName = [cd '\' cutFileName];
end

%% check that filename doesn't already exists
% extInd = strfind(cutFileName,'.');
% counter = 2;
% while exist(cutFileName,'file') == 2
%     warning([cutFileName ' is exisiting file. File name changed to ' cutFileName(1:extInd-1) '_' num2str(counter) '.cut']);
%     cutFileName = [cutFileName(1:extInd-1) '_' num2str(counter) '.cut'];
%     counter = counter + 1;  
% end
if exist(cutFileName,'file') == 2
    [pathCut,nameCut,~] = fileparts(cutFileName);
    if ~isdir([pathCut '\old'])
        mkdir([pathCut '\old']);
    end
    movefile(cutFileName,[pathCut '\old\' nameCut '.cut']);
    warning([cutFileName ' is exisiting file. Old file moved to ' pathCut '\old\' nameCut '.cut']);
end

%% write file
%combine header and clusterIDs for writing
cutFile = [headerStr(:);cellstr(num2str(clustIDs))];
fid = fopen(cutFileName,'w'); %open file for writing
headerLength = length(headerStr);
for i = 1:length(cutFile)
    %header lines with 'min' or 'max' values for centres - need to add
    %white space
    if ~isempty(strfind(cutFile{i},'min')) || ~isempty(strfind(cutFile{i},'max'))
        fprintf(fid,'              %s\r\n',cutFile{i});
    %cluster IDs    
    elseif i > length(headerStr)
        if ~mod(i-headerLength,25) %there is a carriage return every 25 entries
            fprintf(fid,'%s\r\n',cutFile{i});
        else
            fprintf(fid,'%s ',cutFile{i});
        end
    %any other header line    
    else
        fprintf(fid,'%s\r\n',cutFile{i});
    end  
end
fclose(fid);

end

