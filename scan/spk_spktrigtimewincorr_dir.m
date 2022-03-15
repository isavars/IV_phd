function [rtn,binCoordsOut,posMap,shuf] = spk_spktrigtimewincorr_dir(dir, spk_timeA, spk_timeB, T, pos_sample_rate, binsize, varargin)
% Spike-triggered time-windowed spatial correlogram for directional correlates (auto- or cross-).
%
% v2.0 (Aug 2015) incorporates a 'true' spike-time shuffle control (the pos is truely shifted as for the spike times, rather than being a mean displacement histogram).
%
%       [rtn,dirBinVector,posMap,shuffledMapValues] = spk_spktrigtimewincorr(dir, spk_timesA, spk_timesB, timeWinDur, pos_sample_rate, binsize, .. 'paramName', paramVal .. )
%
% Input arguments:
% dir                    - dir position vector for whole trial
% spk_timesA, spk_timesB - spike trains, as a list of spike times in seconds
% timeWinDur             -  is window length, in seconds. If this input is a vector, function will
%                           calculate maps corresponding to all window lengths in the vector, and
%                           return arg will be a cell array, corresponding to these window lengths.
% pos_sample_rate, in Hz
% binsize, in degrees
%
% Optional extra arguments (param name/value pairs):
%
%       'downSample', d       - downsample positions (in time) by a factor of 'd', i.e. take the average dir over d position samples. To reduce memory usage.
%       'dwellThresh',   s    - do not use bins with <s seconds of summed dwell (default 0, i.e. use all possible bins)
%       'mapMode', {'rate', 'pos', 'spk'} - show time windowed rate, or spikes, or position (default='rate').

% Parse optional arguments %
opt.dwellThresh = 0;
opt.mapMode = 'rate';
opt.crossCorrShowReverseTime = 0;
opt.posShuffleActive = 1;
opt.posShuffleMinOffset = 20;
opt.posShuffleStep = 2;
opt.normaliseToShuf = 1;
opt.downSample = 5;
for ii=1:2:length(varargin)
    opt.(varargin{ii}) = varargin{ii+1};
end

% If N pos samples is greater than 32768 (int16 limit), then force a downsample %
if (length(dir) / opt.downSample) > 32768   % Take into account that a downsample may already have been specified.
    opt.downSample = opt.downSample*2;      % Assuming that max trial length is 20 minutes, so downsample by *2 shoudl be sufficient.
end

% If requested, downsample position data in time. Note that 'pos_sample_rate' is also adjusted appropriately at this point. %
if opt.downSample~=1
    d = reshape(dir,opt.downSample,[]);
    dir = rad2deg(   circ_mean(   deg2rad( d ) )  );
    dir = dir';
    dir(dir<0) = dir(dir<0) + 360;
    pos_sample_rate = pos_sample_rate / opt.downSample;
end

spk_sampA = ceil(spk_timeA .* pos_sample_rate);   % TW Note spk_sampA (and spk_sampB) not a histogram, but spike times with sec converted to dwell samp no.
spk_sampB = ceil(spk_timeB .* pos_sample_rate);   %
n_samps = length(dir);
n_spksA = length(spk_sampA);  
win_limit = min( max(T)*pos_sample_rate, n_samps);  % win_limit is the length of the time lag window *in one temporal direction*. Units are number of pos samples, NOT time. Also, if T is vector, win_dur is based on largest value.


% if length(spk_sampA)==length(spk_sampB) && all( (spk_sampA-spk_sampB) == 0 )
%     isautocorr=1;   win_dur = win_limit;
% else
   isautocorr=0;   win_dur = (win_limit*2) + 1;   % win_dur is the length of the whole time lag window, i.e. doubled for cross-corrs
% end


tic
% Unwrap the dir data %
dir = int16( circ_rad2ang( unwrap( circ_ang2rad( dir ) ) ) );

% Work out max size for p_disp_x and pre-allocate
% A note on array format: use 'int16' to try and minimise memeory usage: this means storing dir offsets in degrees and time offsets in pos samples - 
% we are well within the int16 limit for both. The only way to further minimise would be to pre-bin the data, i.e. dir offsets in directional bins and time offsets in sec.
p_disp = zeros( ceil(n_spksA*win_dur), 2, 'int16' );   % Second column contains time lag info, for reconstrcuting different time lag map if T is vector.
s_disp = zeros( ceil(n_spksA * length(spk_sampB) * win_dur/n_samps), 2, 'int16' );


if opt.posShuffleActive
    % Make the shuffled pos array, which is used to construct the shuffled spike maps. This is a 2D array, where each column is the entire dir dataset, but wrapped around a 
    % fixed time point, changing column to column. When you look at each row, therefore, you get the dir when the data has been offset by offset:offset:trialLength values.
    % Then, below in the "1:n_spkA" loop, this will be used to assign 'shuffled' directions to spikes in the windowed correlogram.
    shufOffsets = int16(   opt.posShuffleMinOffset*pos_sample_rate  :  opt.posShuffleStep*pos_sample_rate   :    length(dir) - (opt.posShuffleMinOffset*pos_sample_rate)   );
    shufInd = repmat(shufOffsets,length(dir),1);                % First, each column consists of repeats of the inital offset value.
    shufInd = shufInd + repmat(int16((1:length(dir))'), 1, length(shufOffsets));  % Then add 1:trialLengthInSamps to every column (i.e. starting at initial offset, go to offset + trialLength)
    shufInd(shufInd > length(dir)) = shufInd(shufInd > length(dir)) - length(dir);  % Then subtract trialLengthInSamps from everything greater than trial length, 'wrapping' round data when it goes past the end of the trial.
    shuffledDir = dir(shufInd);
    % Pre-allocate shuffled data arrays. The 'times' arrays are to detect when the windows go over the trial wrapping boundary - the units are pos samp N, rather than sec. %
    p_disp_shuf = zeros(size(p_disp,1), size(shuffledDir,2), 'int16');
    s_disp_shuf = zeros(size(s_disp,1), size(shuffledDir,2), 'int16');  % For the spikes, we will actually produce length(shufOffsets) differernt maps. The divsion of these by the shuf pos map will give us the shuffled rate map population.
%     [p_disp_shuf, p_disp_shuf_times] = deal(zeros(size(p_disp,1), size(shuffledDir,2), 'int16'));
%     [s_disp_shuf, s_disp_shuf_times] = deal(zeros(size(s_disp,1), size(shuffledDir,2), 'int16'));
end



% Construct windowed displacement data %
filled_pvals = 0; filled_svals = 0;  % These three are counter variables for indexing into p_disp and s_disp as they are filled in. We need a separate one for p_shuf: as the shuf pos samples 
pos_win_lengths = nan(1,n_spksA);
for ii = 1:n_spksA;
%     if ii==500;
%         disp('dbstop');
%     end

    ti = spk_sampA(ii);
    if isautocorr;  winStart=ti+1;  else  winStart=ti-win_limit-1;   end
    
    %%% Pos samples %%%
    window_pos = (  max(1,winStart) : min(ti+win_limit, n_samps))';    % 'window_pos' is an index into the vector of pos samples, from the time of the current spkA at 'ti' to the end of the time lag (or the end of the trial). How to fix (for spikes as well)? 
    WL = length(window_pos);    pos_win_lengths(ii) = WL;
    p_disp(filled_pvals+1:filled_pvals+WL, :) = [dir(window_pos)-dir(ti)  int16(window_pos-ti)];  %%% I found the periodic dwell pattern bug! When UNWRAPing here, if the pos window starts *before* pos(ti), pos(ti) won't necessarily be 0, it could be 360, 720 etc... The solution is to unwrap the data only once, before this loop (at a cost of having to use an int16 array).
    % pos offsets for shuffled positions %
    if opt.posShuffleActive
        pDispShufTemp = shuffledDir(window_pos,:)   -    repmat(   shuffledDir(ti,:), length(window_pos), 1 );
        pDispShufTime = shufInd(window_pos,:)   -    repmat(   shufInd(ti,:), length(window_pos), 1 );
        pDispShufTemp( abs(pDispShufTime)>WL ) = 32768;
        p_disp_shuf(filled_pvals+1:filled_pvals+WL, :) = pDispShufTemp;
%         p_disp_shuf(filled_pvals+1:filled_pvals+WL, :) = shuffledDir(window_pos,:)   -    repmat(   shuffledDir(ti,:), length(window_pos), 1 );
%         p_disp_shuf_times(filled_pvals+1:filled_pvals+WL, :) = shufInd(window_pos,:)   -    repmat(   shufInd(ti,:), length(window_pos), 1 );
    end
    filled_pvals = filled_pvals + WL;  % Bump pos counter. Note even the last interation is important, used to crop the pre-allocated matrix, below.

    
    %%% Spikes %%%
    window_spk = spk_sampB( spk_sampB>=winStart & spk_sampB <= ti+win_limit );   % Get the spikes (B) in the time window (units=pos samp N)
    WL = length(window_spk);
    if WL==0;   continue;   end;   % If there are no spikes in window
    s_disp(filled_svals+1:filled_svals+WL, :) = [    dir(window_spk)-dir(ti)      int16(window_spk-ti)  ];    % To get the spike positions, index into the existing unwrapped dir for the window_spk
    % Spike offsets for shuffled positions %
    if opt.posShuffleActive
        sDispShufTemp = shuffledDir(window_spk,:)   -    repmat(   shuffledDir(ti,:), length(window_spk), 1 );
        sDispShufTime = shufInd(window_spk,:)   -    repmat(   shufInd(ti,:), length(window_spk), 1 );
        sDispShufTemp( abs(sDispShufTime) > length(window_pos) ) = 32768;
        s_disp_shuf(filled_svals+1:filled_svals+WL, :) = sDispShufTemp;
    end
    filled_svals = filled_svals + WL;  % Bump spike counter. Note even the last interation is important, used to crop the pre-allocated matrix, below.
    
end



% strip out the pre-allocated values
p_disp = p_disp( 1:filled_pvals , : );     
s_disp = s_disp( 1:filled_svals , : );
if opt.posShuffleActive
    p_disp_shuf = p_disp_shuf( 1:filled_pvals , : );
    s_disp_shuf = s_disp_shuf( 1:filled_svals , : );
end

% Set 'backwards-in-time' values to be negative .
if  isautocorr
    % add each disp and its negative to count each pair in each order %
    p_disp = [p_disp; -p_disp]; 
    s_disp = [s_disp; -s_disp];
    p_disp(:,2) = abs(p_disp(:,2));
    s_disp(:,2) = abs(s_disp(:,2)); % Time lag data should stay absolute in this case
elseif opt.crossCorrShowReverseTime
    p_disp( p_disp(:,2)<0,  1   ) = -p_disp( p_disp(:,2)<0,  1   );
    s_disp( s_disp(:,2)<0,  1   ) = -s_disp( s_disp(:,2)<0,  1   );
    p_disp_shuf( p_disp_shuf(:,2)<0,  1   ) = -p_disp_shuf( p_disp_shuf(:,2)<0,  1   );
elseif ~opt.crossCorrShowReverseTime
    % When forward and reverse time are collapsed %
    p_disp( :,  2  ) = abs(  p_disp( :,  2  )  );
    s_disp( :,  2  ) = abs(  s_disp( :,  2  )  );
%     p_disp_shuf( :,  2  ) = abs(  p_disp_shuf( :,  2  )  );
end

% Make maps (histograms of pos and spiking) %
if isautocorr || ~opt.crossCorrShowReverseTime;    
    timeLagLims = [zeros(size(T)); T] * pos_sample_rate;
else
    timeLagLims = [[fliplr(-T); zeros(size(T))], [zeros(size(T)); T]] * pos_sample_rate;     % If requested also create the 'backwards-in-time' rate maps.
end;
binVect = floor(min(p_disp(:,1))/binsize)*binsize : binsize : ceil(max(p_disp(:,1))/binsize)*binsize;  % One bin vector for all time lags, based on absolute max and min displacement. Makes aligning of time lag maps (in calling functions) much easier.


% Pre-assign maps (these will be multi-row for when T is a vector, one row per time-lag) %
[posMap,spkMap] = deal(nan(length(T),length(binVect)));
s = {'binMean','binSD','bin95'};
for ii=1:length(s);   shuf.(s{ii}) = nan(length(T),length(binVect));   end  % Assign the shuf values anyway, so that 
s = {'timeLagMean','timeLagSD','timeLag95'};
for ii=1:length(s);   shuf.(s{ii}) = nan(length(T),1);   end

% Loop through Time Lags %
for ii=1:size(timeLagLims,2)
    % Filter displacement data so as to get only the appropriate time lags
    p_disp_forTimeLag = p_disp( p_disp(:,2)>timeLagLims(1,ii) & p_disp(:,2)<=timeLagLims(2,ii), 1);
    s_disp_forTimeLag = s_disp( s_disp(:,2)>timeLagLims(1,ii) & s_disp(:,2)<=timeLagLims(2,ii), 1);
    % Make Rate Map %
    posMap(ii,:) = transpose( histc(p_disp_forTimeLag,binVect) );
    spkMap(ii,:) = transpose( histc(s_disp_forTimeLag,binVect) );
    
    if opt.posShuffleActive
        % Make shuffled rate maps from spike counts %
        posMapShuf = histc( p_disp_shuf(          p_disp(:,2)>timeLagLims(1,ii) & p_disp(:,2)<=timeLagLims(2,ii),:)          , binVect , 1);
        spkMapShuf = histc( s_disp_shuf(          s_disp(:,2)>timeLagLims(1,ii) & s_disp(:,2)<=timeLagLims(2,ii),:)          , binVect , 1);  % For the shuffled spike map, the time lag index is the same as for normal data.

        rateMapShuf = (spkMapShuf ./ posMapShuf) .* pos_sample_rate;  % Convert to Hz as we are making the rate map.
        rateMapShuf( posMapShuf < opt.dwellThresh*pos_sample_rate ) = nan;
        
        % Get the various measures of dispersion from the shuffled maps %
        rateMapShuf = rateMapShuf';  % Switch in to the same orientaion as the multi-lag results arrays (column = dir offset).
        shuf.binMean(ii,:) = nanmean(rateMapShuf,1);
        shuf.binSD(ii,:) = nanstd(rateMapShuf,1);
        shuf.bin95(ii,:) = prctile(rateMapShuf,95,1);
        shuf.timeLagMean(ii,1) = nanmean(rateMapShuf(:));
        shuf.timeLagSD(ii,1) = nanstd(rateMapShuf(:));
        shuf.timeLag95(ii,1) = prctile(rateMapShuf(:),95);
    end
    
    if size(spkMap,1)==1;   spkMap=spkMap';   end  % This happens when n spike in map = 1.
    
end
% Filter for points of low dwell time %
spkMap(posMap < opt.dwellThresh*pos_sample_rate) = nan;
posMap(posMap < opt.dwellThresh*pos_sample_rate) = nan;
rateMap = (spkMap./posMap) .* pos_sample_rate;   % Convert to Hz as we are mkaing the rate map.

toc

if opt.posShuffleActive && opt.normaliseToShuf
    rateMap = (rateMap - shuf.binMean)  ./  shuf.binSD;
%     rateMap(rateMap<3) = 0;
end

rtn = rateMap;

% Rate map or spk map %
% if strcmp(opt.mapMode,'rate')
%     rtn = spkMap./posMap;
% elseif strcmp(opt.mapMode,'pos')
%     rtn = posMap;
% elseif strcmp(opt.mapMode,'pos')
%     rtn = spkMap;
% end

binCoordsOut = binVect;

clear p_disp s_disp full_p_disp full_s_disp
if opt.posShuffleActive
    clear s_disp_shuf p_disp_shuf posShufInd shufInd shuffledDir
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [sm_map]=smoothPlaceMap(map,k,unvis)
% % Smooth a place map. map=map, k=boxcar kernel, unvis=index of unvisited bins.
% if max(size(k))==1;  k=ones(k);  end   % Expand single parameter to flat k-by-k square
% map(unvis)=0;
% visTemplate=ones(size(map));
% visTemplate(unvis)=0;
% filtMap=imfilter(map,k);
% filtVis=imfilter(visTemplate,k);
% sm_map=filtMap./filtVis;
% sm_map(unvis)=nan;

















