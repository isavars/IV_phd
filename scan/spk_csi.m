function [res, res2] = spk_csi(cellData)
%   [csi, csi_wilson] = spk_csi(cellData);
%
% csi & csi_wilson: whats the difference?
% csi - based on proportion of spikes that have an appropriate inter-spike interval.
% csi_wilson - based on proportion of total number of spikes.
%
% csi makes more sense, csi_wilson is how it was in the code he gave to us.

if (any(length(cellData.st) == [0 1]))  
    res = NaN;  res2 = NaN;   return
end

min_interval = 0.003;
max_interval = 0.015;

% Parse input %
if cellData.wf_mode==2
    amps=double(cellData.wf_amps);
    [dum,maxCh] = max(mean(amps,1));                   % Find channel with highest mean amplitude. (2nd argout gives index)
    amps = amps(:,maxCh);                              % amp is amplitude for each spike in the channel with the highest mean amp.
elseif cellData.wf_mode==1;
    wf = double(cellData.wf_all);
    amps = squeeze( max(wf,[],2) - min(wf,[],2) )';     % Find amplitudes for every spike on every channel.
    [dum,maxch] = max(mean(amps,2));                    % Find channel with highest mean amplitude. (2nd argout gives index)
    amps = amps(maxch,:)';                              % amp is amplitude for each spike in the channel with the highest mean amp.
else
    disp('Can''t do CSI - no amplitude data!');
    res = NaN;  res2 = NaN;   return
end
times=cellData.st;
    
% Put times and WF into correct format for original wilson code %
if length(times) == length(amps)
    for i=1:length(times)
        point(i).h = amps(i);
        point(i).t = times(i);
    end
else
    error('Different numbers of spike times and spike amplitude');
end


% These are for testing: comparison to regular version %
sub = 0;
tl = 0;
tr = 0;
bl = 0;
br = 0;

% Complex spike index %
csi = 0;
ncsi = 0;
for i=2:(length(point) - 1)
    % Get interval for this spike %
    interval1 = point(i).t - point(i-1).t;
    interval2 = point(i).t - point(i+1).t;
    
    % Get smallest interval %
    if interval1 < -interval2
        ht = point(i).h - point(i-1).h;
        minint = interval1;
    else
       ht = point(i).h - point(i+1).h;
        minint = interval2;
    end
    
    % Compute CSI %
    if abs(minint) <= max_interval
        if abs(minint) < min_interval
            csi = csi - 1;
            sub = sub + 1;
        elseif (ht <= 0) && (minint > 0)
            csi = csi + 1;
            br = br+1;
        elseif (ht > 0) && (minint < 0)
            csi = csi + 1;
            tl = tl+1;
        elseif (ht <= 0) && (minint < 0)
            csi = csi - 1;
            bl = bl+1;
        elseif (ht > 0) && (minint > 0)
            csi = csi - 1;
            tr = tr+1;
        end
        ncsi = ncsi + 1;
    end
    
end
    
res = (csi / ncsi) * 100;                   % this one seems to make more sense
res2 = (csi / (length(point) - 2)) * 100;   % This is the one in wilsons code
    
