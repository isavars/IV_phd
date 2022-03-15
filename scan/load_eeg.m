function [rtn] = load_eeg(trialName)

% Load EEG data from DACQ file
%
%       [rtn] = load_eeg('trialName.eeg');
%
% rtn.eeg           - voltage samples (int8)
% rtn.sample_rate   - Sample rate.
%
% Note that file extension should be included in input argument (as these
% change depending on Dacq SW Version).
%
% %% TODO (05/09/08): This only reads .eeg files (250Hz), not .egf files
% (4800Hz). Should do both. Watch for 1) different header formats, 
% 2) 2 bytes/sample.

% Preallocate dummy output in case of data read faliure.
rtn.eeg = [];        % If we have one of these, just quit now.
rtn.sample_rate = [];

% Read in data %
fid = fopen(trialName,'r','ieee-be');  % 'ieee-be' is machine format, 'big endian'.
if fid==-1
    fprintf(1, 'SCAN: Problem loading EEG. Can''t find file %s \n', trialName);
    return
end

hdr = fread(fid,400,'int8');
ds  = strfind(hdr','data_start') + 10; % data start marker 
% de  = strfind(hdr','data_end') - 3; % data end marker

% more convenient format
frewind(fid);
tempHeader = textscan(fid,'%s %[^\n]',11);
tempHeader = horzcat(tempHeader{:});

sRateInd = strcmp('sample_rate',tempHeader(:,1));
if sum( sRateInd )==0   % Sometimes there are weird 'empty' EEG files, which don't even have a full header. If we have one of these, just quit now.
    fprintf(1, 'SCAN: Problem loading EEG. No EEG data contained in file %s \n', trialName);
    return
end
rtn.sample_rate = sscanf(tempHeader{sRateInd,2},'%d');

bytesInd = strcmp('bytes_per_sample',tempHeader(:,1));
bytesPerSamp = sscanf(tempHeader{bytesInd,2},'%d');

fseek(fid,ds,'bof');  
if bytesPerSamp==1
    nSamplesInd = strcmp('num_EEG_samples',tempHeader(:,1));
    nSamples    = sscanf(tempHeader{nSamplesInd,2},'%d');
    rtn.eeg = fread(fid,nSamples,'int8');
%     rtn.eeg = int8(binData(ds:de));   % simple
elseif bytesPerSamp==2
    nSamplesInd = strcmp('num_EGF_samples',tempHeader(:,1));
    nSamples    = sscanf(tempHeader{nSamplesInd,2},'%d');
    %grab actual voltage data 
    rtn.eeg = fread(fid,nSamples,'int16'); %re-read as int16
    
%     rtn.eeg = int16(data(1:nSamples));
    
end
fclose(fid);

end





% 
% % Read Header: old-fashioned style %
% headerText = fread(fid,800,'int8');
% headerText = char(abs(headerText))';
% 
% sr_index = (findstr('bytes_per_sample ',headerText))+16;     % Get number of bytes per sample (is it a .eeg or a .egf file?)
% bytesPerSamp = str2double( strtok(headerText(sr_index:end)) );
% 
% sr_index = (findstr('sample_rate ',headerText))+11;     % Get sample rate
% sr = strtok(headerText(sr_index:end));
% rtn.sample_rate = str2double(sr);
% 
% ds = (findstr('data_start',headerText))+10;        % Look for data start marker
% 
% % Get data %
% if bytesPerSamp==1
%     m = (findstr('num_EEG_samples ',headerText))+16;       % Look for number of samples marker
% %     nsamp = str2double(headerText(m:(m+9)));               % get number of samples
%     nsamp = sscanf(headerText(m:(m+9)),'%d');                 % LM bug fix: when making joined files the lack of the whitespace causes issues
%     % bug with joined files; LMfix
%     %     fseek(fid,(ds),'bof');                               % Set file position indicator at data start 
%     frewind(fid)
%     data = fread(fid,'int8');               % Read data (header and voltage samples). 'int8' converts each byte into a integer -126:126.
%     rtn.eeg = int8(data(ds:ds+nsamp-1)); %
% %     rtn.eeg = int8(data(1:nsamp));
%   
% else
%     m = (findstr('num_EGF_samples ',headerText))+16; % Look for number of samples marker
%     nsamp = sscanf(headerText(m:(m+9)),'%d');  % LM bug fix: when making joined files the lack of the whitespace causes issues
% %     nsamp = str2double(headerText(m:(m+9)));               % get number of samples
% %     fseek(fid,ds,'bof');                               % Set file position indicator at data start 
%     frewind(fid)
%     data = fread(fid,'int16');     
%     rtn.eeg = int16(data(ds:ds+nsamp-1)); %
% %     rtn.eeg = int16(data(1:nsamp));
%     
% end
% 
% 
% fclose(fid);