function [rtn] = load_pos(trialName, nLight)

% Load in positions from a .pos file
%
%       rtn = load_pos('trialName', nLight);
%
% rtn.led_pos     - positions, in format [nSamp, nLight, x:y]
% rtn.led_pix     - LED size in pixels, format [nSamp, nLight]                  
% rtn.sample_rate
% rtn.num_pos_samples
% rtn.pos_format
% rtn.header      - Header in {n,2} cellstr format, for use with mTint KEY_VALUE.
%
% NOTE: the fieldnames and format of the first two are for compatibility
% with POSTPROCESS_POS_DATA (mTint).

%%% Read data from File %%%
fid = fopen([trialName '.pos'],'r','ieee-be');
if fid==-1
    ME = MException('scan:load_pos:posFileNotFound', ['Could not open pos file ' fileName '.pos']);
    throw(ME);
end
% Read Header: old-fashioned style %
headerText = fread(fid,800,'int8');
headerText = char(abs(headerText))';
headerFields = {'sample_rate', 'num_pos_samples', 'pos_format', 'pixels_per_metre', ...
                'window_max_x', 'window_min_x', 'window_max_y', 'window_min_y'};
for ii=1:length(headerFields);
    ind = (findstr(headerFields{ii},headerText))+length(headerFields{ii});
    if strcmp(headerFields{ii},'pos_format')
        header.(headerFields{ii}) = strtok(headerText(ind:end));
    else
        header.(headerFields{ii}) = str2double(strtok(headerText(ind:end)));
    end
end
ds = (findstr('data_start',headerText))+10;                 % Look for data start marker
% Read data %
fseek(fid,(ds-1),'bof');                                    % Set file position indicator at data start 
data = fread(fid, header.num_pos_samples*10, 'uint16');
fclose(fid);

%%% Reshape into correct output format %%%
data = reshape(data, [10 header.num_pos_samples])';
data = data(:,3:10);
data = reshape(data, [header.num_pos_samples, 2, 4]); % Now in format [nSamp, x:y, nLight]
% Separate numpix, if existing, switch format to [nSamp, nLight, x:y] %
if strcmp(header.pos_format, 't,x1,y1,x2,y2,numpix1,numpix2') && nLight<=2
    led_pos = repmat(NaN, [header.num_pos_samples, nLight, 2]);
    led_pix = repmat(NaN, [header.num_pos_samples, nLight]);
    for ii=1:nLight
        for jj=1:2
            led_pos(:,ii,jj) = data(:,jj,ii);
        end 
        led_pix(:,ii) = data(:,ii,3); % numpix always seems to start at 3rd light (5th entry)
    end
else
    led_pos = repmat(NaN, [header.num_pos_samples, nLight, 2]);
    led_pix = [];
    for ii=1:nLight
        for jj=1:2
            led_pos(:,ii,jj) = data(:,jj,ii);
        end
    end
end
led_pos(led_pos==1023) = NaN;   % mTint functions are 
led_pix(led_pos==1023) = NaN;   % expecting this.

% Read header: mTint format (needed for compatibility) %
[key value]=textread([trialName '.pos'],'%s %[^\n]',27); % CAUTION! N Lines of header hard-coded.
textHeader=[cat(1,key) cat(1,value)];

% Output %
rtn = header;
rtn.header = textHeader; % For mTint function compatibility.
rtn.led_pos = led_pos;  
rtn.led_pix = led_pix;
