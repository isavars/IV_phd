
function [output] = load_spikes(filename)

% Reads a tetrode file (ie .1 .2 etc.) and returns the spike waveforms and times.
% 
%         `spikes.times      = read_spikes('filename');
%                .waves
%                .trialDur
%         
% - spikes.times is a vector containing spike times, first to last. Units are seconds
% - spikes.wave is a 3-D matrix, with format [spike number(first-last), voltage samples(1-50), channels(1-4)]
% - spikes.trialDur is the trial duration, in seconds.


fid = fopen(filename,'r','ieee-be');        % Open tet file. 'ieee-be' string is machine format, which needs to be 'big endian' to read times correctly. Default format doesnt work on my PC.
if fid==-1
    error(['SCAN: Error in LOAD_SPIKES. Can''t find file ', filename])
end
data = fread(fid,'int8');               % Read data (header and voltage samples). 'int8' converts each byte into a integer -126:126.

%%%%%

hdr = char(abs(data(1:400)))';              % Translate header into text. Use 'abs' to avoid error message on negative integers after data starts.

m = (findstr('num_spikes ',hdr))+11;        % Look for number of spikes marker
% nspk = str2double(hdr(m:(m+9)));               % get number of spikes
nspk = sscanf(hdr(m:(m+9)),'%d');                 % LM bug fix: when making joined files the lack of the whitespace causes issues

trial_dur_index = (findstr('duration',hdr))+8;      % Get trial duration.
trial_dur = strtok(hdr(trial_dur_index:end));
trial_dur = str2double(trial_dur);

timebase_index = (findstr('timebase ',hdr))+8;     % Get sample rate (DSP timer Ticks/s)
timebase = strtok(hdr(timebase_index:end));
timebase = str2double(timebase);

ds = (findstr('data_start',hdr))+10;        % Look for data start marker

%%%%%%
                                                
data = data(ds:(ds+(nspk*216)-1));           % Take voltage sample values. (Time stamps still included)
data = reshape(data,54,4,nspk);               % Reshape into (sample,channel,spike)
data = data(5:end,:,:);                       % Cut off time stamp samples
data = shiftdim(data,2);                      % Make 'spikes' first(row) dimension of output matrix.
data = int8(data);                            % Convert to int8 (-127 - 128) format.


%%%%%%


fseek(fid,(ds-1),'bof');                        % Set file position indicator at data start (for reading time stamps, see below)
times = fread(fid,nspk,'int32',212);            % read time stamps. 'nspk'=read this many samples, 'int32'=read 32bit integer, '212'=skip 212 bytes between reading values.


%%%%%%

output.times = times ./ timebase;
output.waves = data;
output.timeBase = timebase;
output.trialDur = trial_dur;

fclose(fid);



