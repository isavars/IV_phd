function make_tint_file_gain_troubleshoot(tetrode_mat,tetFileName,tetrode_number,trialInfo) %,DACQinfo -removed
%% Generate TINT file stucture
% Timestamp = 4 byte timestamp (big-endian)
% Channel samples = 50 8-bit samples
% Total = 54 bytes per channel

% % Tetrode filenames 1-4 or 5-8 (Left/Right plug)
% if strcmp('Right',trialInfo.plug)
%     % DACQ tetrodes 5 -8
%     tetrode_number = tetrode_number + 4; 
% else
%     %DACQ tetrodes 1 - 4
% end

%% Convert timestamp to 4 bytes and big-endian
timebase = 96000;
time = tetrode_mat(:,1)'; % convert to row vector
time = time./1e6; % Convert to seconds

time = int32(time.* timebase); % Convert to int32 and multiply by timebase

% Convert timestamp to big-endian if needed (Windows PC = little-endian)
[~,~,endian] = computer;
if strcmp(endian,'L')
    time = swapbytes(time);
end

time = reshape( typecast(time,'int8'),4,[]);

%% Convert signal to 50 8-bit samples
voltages = tetrode_mat(:,2:end)';

scaleMax = 100; % can hardcode because its always the same for me  but check if this value makes sene - DACQinfo.scaleMax(tetrode_number,:);
DACQ_sample = nan(size(voltages));
% Adapt conversion of each channel depending on gain used in DACQ
for i = 1:4
    DACQ_sample(:,i:4:end) = ceil((voltages(:,i:4:end)/scaleMax.*128)); %removed index (i) from scalemax because it's not needed 
end
      
% DACQ_sample = ceil((voltages/scaleMax) .* 128);
DACQ_sample(DACQ_sample == -128) = -127;

% Convert to int8 datatype
DACQ_sample = int8(DACQ_sample);

%% Combine timestamps and voltage data 

final_vec = int8(zeros(54,size(DACQ_sample,2)));
final_vec(1:4,:) = time;
final_vec(5:end,:) = DACQ_sample; % append voltage data

final_vec = reshape(final_vec,216,[]); % 1 row = 1 spike, all channels

%% Write to file 

fid = fopen([trialInfo.writeDir '\' tetFileName],'a+');
fprintf(fid,'%s','data_start'); %add data marker
fwrite(fid,final_vec(:),'int8',0,'ieee-be');
% fwrite(fid,final_vec,'int8'); % dtype = int8
fprintf(fid,'\r\n%s\r\n','data_end'); % need at end of Tint file
fclose(fid);

fprintf(strcat("Finished generating tet file for Tetrode ",num2str(tetrode_number),'\n'));







