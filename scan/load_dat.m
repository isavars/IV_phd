
function [rtn] = load_dat(trialName)

% Read position and EEG data from .dat file.
%
%           rtn.x = load_dat('trialName');
%              .y
%              .dir
%              .speed     
%              .sample_rate
%              .ppm
%
% Output format for each (data) field of rtn is (first_sample:last_sample).


% Check if data files exist %   
% dat_namer_struct = dir([trialName, '_*.dat']);	
% if isempty(dat_namer_struct)
%     error(['Could''t find a .dat file for trial ', trialName]);
% elseif length(dat_namer_struct) == 1
%     dat_namer = dat_namer_struct(1).name;
% elseif length(dat_namer_struct) > 1
%     dat_namer = dat_namer_struct(1).name;
%     fprintf(1,'\n%s', ['(More than one dat file for ', trialName, '. No problem, but will only use ', dat_namer, ').']);
% end

dat_namer = trialName;

%%%%% Read Header %%%%% 
% NOTE: I know this seems like a long-winded way to read a header, but I could not get fucking
% textread, strread, or fscanf to work. textread and strread seem to have a genuine bug (give a error
% for any template not on first line), and by then, I was too annoyed to try and understand why
% fscanf wasnt working.

fid = fopen(dat_namer);
hdr = fread(fid,300);
hdr = char(hdr');
fclose(fid);

temp_ind = findstr(hdr,'converted at (')+14;
ppcm = str2double(strtok(hdr(temp_ind:end),44));    % strtok to comma
temp_ind = findstr(hdr,'start')+6;
start_time = str2double(strtok(hdr(temp_ind:end)));
temp_ind = findstr(hdr,'end')+3;
end_time = str2double(strtok(hdr(temp_ind:end)));
temp_ind = findstr(hdr,'rows')+4;
rows = str2double(strtok(hdr(temp_ind:end)));

% Get positions %
data = dlmread(dat_namer,'\t',[5,1,(rows+4),6]);
positions = data(:,1:4);
    
% Check sample rate, if necessary, downsample (for eeg sample rate) %
samp_rate = rows/(end_time-start_time);
pos_samp_rate = samp_rate;
if 0%samp_rate > 50
    % get true sample rate %
    fid = fopen([trialName, '.pos']);
    hdr = char(fread(fid, 400)');
    fclose(fid);
    temp_ind = findstr(hdr,'sample_rate')+12;
    pos_samp_rate = str2double(strtok(hdr(temp_ind:end)));
    % downsample index %
    ds_step = samp_rate / pos_samp_rate;
    ds_start = ceil(ds_step/2);
    ds_end = size(positions,1) - floor(ds_step/2);
    ds_ind = ceil(ds_start:ds_step:ds_end);
    % downsample %
    positions = positions(ds_ind,:);           
elseif samp_rate < 46.875
    disp(['Problem with trial ', trialName, '. Dat file sample rate too low to read positions.']);
end
       
% Convert X+Y from cm to pixels %
positions(:,1) = round(positions(:,1).*ppcm);
positions(:,2) = round(positions(:,2).*ppcm);
 
%%% Assign to output structure %%%
rtn.sample_rate = pos_samp_rate;
rtn.ppm = ppcm*100;
rtn.x = positions(:,1);
rtn.y = positions(:,2);
rtn.dir = positions(:,3);
rtn.speed = single(positions(:,4));
rtn.eeg = data(:,5);
rtn.phase = data(:,6);
% if samp_rate >= 234
%     rtn.eeg = data(:,5);
%     rtn.phase = data(:,6);
%     rtn.eegSampRate = samp_rate;
% else
%     rtn.eeg = [];
%     rtn.phase = [];
%     rtn.eegSampRate = [];
% end
    






