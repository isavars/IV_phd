function [rtn,binCoordsOut,shuf] = spk_spktrigtimewincorr(x, y, spk_timeA, spk_timeB, T, pos_sample_rate, binsize, varargin)
% Spike-triggered time-windowed spatial correlogram.
% Incorporates a spike-time shuffle control (spike train B is offet relative to position, a user-defined number of times).
%
%       [rateMap,distBinVector] = spk_spktrigtimewincorr_displace(X, Y, spk_timesA, spk_timesB, timeWinDur, pos_sample_rate, binsize);
%       [rateMap,distBinVector] = spk_spktrigtimewincorr_displace(DIR, [], spk_timesA, spk_timesB, timeWinDur, pos_sample_rate, binsize);
%
%       [rateMap,distBinVector,shufData] = spk_spktrigtimewincorr_displace(X, Y, spk_timesA, spk_timesB, timeWinDur, pos_sample_rate, binsize, .. 'paramName', paramVal .. );
%       [rateMap,distBinVector,shufData] = spk_spktrigtimewincorr_displace(X, Y, spk_timesA, spk_timesB, timeWinDur, pos_sample_rate, binsize, paramStruct );
%
%       
%
% FIXED INPUTS:
%
% X, Y                    - Position data. X and Y position vectors *OR* DIR vector, plus [] as a placeholder.
%                           !!!VERY IMPORTANT!!! Note that, to save memory, x and y are discretised inside the function, so large units
%                           (e.g. cm) will be rounded, with accompanying loss of accuracy. !RECOMMENDED TO USE PIXELS OR MM AS POS UNITS!
%                           (for DIR data use degrees, not radians - this is assumed inside the function).
% spk_timesA, spk_timesB  - spike trains, as a list of spike times in seconds
% timeWinDur              - Window length, in seconds.
% pos_sample_rate         - For the x and y data, in Hz
% binsize                 - bin size for rate maps, in pixels
%
% OUTPUTS:
%
% 'rateMap'             - This is the 'correlogram', actually, a time-windowed rate map of spk B firing relative to spk B. 
%                         Units are Hz or Z-scores of shuffled control, depending on options specified.
% 'distBinVector'       - Vector of spatial displacment bins, corresponing to the rates in 'rateMap'. The units are the same as input x or y.
% 'shufData.binMean'    - Mean of the shuffled control rate maps. Bins correspond to those in 'rateMap' and 'distBinVector'.
% 'shufData.binSD'      - S.D. of the shuffled control rate maps.
% 'shufData.maps'       - The actual shuffled maps, an ( length(binVect),length(binVect),nPosShuffles ) array, where each 2D array (:,:,n) is one shuffled map.
%
% OPTIONAL PARAMETERS:
% (These can be specified as param name/value pairs in a comma-separated list, or, if you prefer, you can pass a single structure, where
% the field has the name of the argumant and the value of the field is the argument).
%
%  'temporalDownSamp', d        - temporal down-sample positions (in time) by a factor of 'd', i.e. take the average displacement over d position samples. To reduce memory usage.
%  'dwellThresh',   s           - exclude bins with <s seconds of summed dwell from the rateMap output (default 2sec)
%  'boxcar', N                  - N-bin long boxcar filter width for rate map smoothing. If set to [], no smoothing done.
%  'posShuffleActive', 1 or 0   - Activate position shuffle control (this is memory intensive, so basically there is the option to turn it off if you are not using it).
%  'normaliseToShuf', 1 or 0    - If 1, the 'rateMap' output is converted to Z-scores of the shuffled population (i.e. S.D.s above the mean).
%  'posShuffleMinOffset, S      - Minimum position offset for shuffled control. Default 20s, I don't see any reason to change this.
%  'nPosShuffles',   N          - Number of pos shuffles in control. Be careful of memory limits here, as all of the big variables in the function scale linerarly 
%                                 in size with this parameter. 100 is OK to get a rough but correct answer, I think, 500 for super-reliable convergence. 200 for normal use?
%  'offsetTimeWindow', 1 or 0   - Run the function with a time window that is (spkA+T1):(spkA+T2). (Whereas normal use is spkA:(spkA+T)). Need to set T input
%                                 to a 1x2 column vector, [time1 time2], in seconds.
%  'showAs1D', 1 or 0           - Collapse x and y to 1D radial displacement.
%  'longestDist', 1 or 0        - Collapse x and y to 1D displacement, 'longest distance', i.e. acutal summed movement along path. WARNING: don't use this mode with pix as units,
%                                 you will get overflow of the int16 datatype. Convert to mm first, instead.


% Some version history:
% v2 - incorporating DMs code for accurately pre-allocating arrays, re-instating multiple maps from one call, and possbibly 2D maps.
%  a - DM code for making code more efficent is done, also some bsxfun in my code.
%  b - change map loop code for multiple time lags. I did some preparatory work for 2D maps (harmonising array dimensions), but haven't done the tricky bit, changing histc to accumarray for map building.
%  c - Change the time window specification so that you can run the window (spkA+time1 : spkA+time2). (Previously limited to (spkA : spkA+time1))
%  d - 2D map output
%
% Still to be incorporated into this function (03/Mar/2016):
% - directional data
% - reinstate the option to get multiple time lags from one run of the function.
%
% UPDATE 05/May/2016
% - dir data done and working.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sampling, smoothing, downsampling %
opt.dwellThresh = 2;
opt.boxcar = [];
opt.temporalDownSamp = 5;
% Shuffle control %
opt.posShuffleActive = 1;
opt.posShuffleMinOffset = 20;
opt.nPosShuffles = 200;
opt.normaliseToShuf = 1;
% Time window options %
opt.offsetTimeWindow = 0;
% Running mode for the function. These three parameters shoudl really be collapsed into one,
% as they are all mutually exclusive. When they are all set to zero, the 'default' running mode
% is the 2D (x,y) verseion of the TW-CC.
opt.showAs1D    = 0;    % Collapse x and y to 1D radial displacement, shortest distance 'as-the-crow-flies'.
opt.longestDist = 0;    % Collapse x and y to 1D, but now take the longest distance along the path (i.e. actual summed movement).
opt.isDirData   = 0;    % Directional version.
if isempty(y);   opt.isDirData = 1;    else   opt.isDirData = 0;   end  % Also autodetect directional version from input arguments, so caller doesn't have to specifiy.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input and various bits of set-up %
if ~isempty(varargin)
    if ischar(varargin{1})
        for ii=1:2:length(varargin)
            opt.(varargin{ii}) = varargin{ii+1};
        end
    elseif isstruct(varargin{1})
        s=varargin{1};  f=fieldnames(s);
        for ii=1:length(f)
            opt.(f{ii}) = s.(f{ii});
        end
    end
end
% The type of the big offset arrays will depend on the function running options %
castFun = @int16;   typeString = 'int16';   % NOTE: actually, this is not true at the moment, INT16 is used for everything, but the code is left in place anyway.
% The size of big arrays will depend on running mode:
% A quick note on why the 'extra' columns (time, nRepShuf) are necessary (see below): time is required in order to be able to calculate the maps
% for mulitple time lags from only one run of the function, nRepShuf is implicit in the order of the 3rd dimension of x_disp_shuf, but needs to 
% be supplied as an explicit index to ACCUMARRAY in order to make 1x1xnRepShuf shuffled maps.
if opt.isDirData  ||  opt.longestDist
    nDataDims = 1;          % Dir data or summed displacement data: x_disp has two columns, dir/disp and time
    nMapDims = 1;           % (for summed displacment, x,y) is collapsed to 1D 'distance' *before* spike loop.
    nColInShufArray = 2;    % x_disp_shuf has two columns, dir/disp and nRepShuf
elseif opt.showAs1D
    nDataDims = 2;          % 2D data collapsed to radial displacement:  x_disp has three columns, x, y and time
    nMapDims = 1;           %
    nColInShufArray = 2;    % x_disp_shuf has two columns, these are initially x and y, during the main loop, but then become radDisp and nRepShuf afterwards, before making maps.
else
    nDataDims = 2;          % 2D data:  x_disp has three columns, x, y and time
    nMapDims = 2;           %
    nColInShufArray = 3;    % x_disp_shuf has three columns, x, y and nRepShuf
end
% Check and amend format of T, depending on opt.offsetTimeWindow mode %
if opt.offsetTimeWindow 
    if numel(T)~=2
        error('If running function in ''offsetTimeWindow'' mode, T should be a 2-element vector ([Time1 Time2])');
    elseif size(T,1) > size(T,2)
        T=T';   % In this case, format of T is [innerWinLimit, outerWinLimit].
    end
end
if ~opt.offsetTimeWindow
    T = [zeros(length(T),1) T(:)];   % if running in 'normal' mode (window runs spkA:spkA+T(ii)), fill in the implicit zero for the length of the inner window 
end                                  % and set format to [0 timeLag1; 0 timeLag 2; 0 timeLag 3; etc].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the data for making AC/CC %
% If N pos samples is greater than 32768 (int16 limit), then force a downsample %
if (length(x) / opt.temporalDownSamp) > 32768   % Take into account that a downsample may already have been specified.
    opt.temporalDownSamp = opt.temporalDownSamp*2;      % Assuming that max trial length is 20 minutes, so downsample by *2 shoudl be sufficient.
end
% If requested, downsample position data in time. Note that 'pos_sample_rate' is also adjusted appropriately at this point. %
if opt.temporalDownSamp~=1
    if opt.isDirData
        x = reshape(double(x), opt.temporalDownSamp, [] );         % Dir pos data comes into function as 'x' input.
        x = rad2deg(   circ_mean( deg2rad(x), [], 1 )  )';
        x(x<0) = x(x<0) + 360;
    else
        x = nanmean( reshape(double(x),opt.temporalDownSamp,[]), 1 )';
        y = nanmean( reshape(double(y),opt.temporalDownSamp,[]), 1 )';
    end
    pos_sample_rate = pos_sample_rate / opt.temporalDownSamp;
end
% Prepare some general variables %
spk_sampA = ceil(spk_timeA .* pos_sample_rate);   % TW Note spk_sampA (and spk_sampB) not a histogram, but spike times with sec converted to dwell samp no.
spk_sampB = ceil(spk_timeB .* pos_sample_rate);   %
n_samps = length(x);
n_spksA = length(spk_sampA);    n_spksB = length(spk_sampB);  
win_limit = min( max(T(:,2))*pos_sample_rate, n_samps);        % win_limit is the length of the time lag window *in one temporal direction*. Units are number of pos samples, NOT time. Also, if T is vector, win_dur is based on largest value.
win_limit_inner = min( min(T(:,1))*pos_sample_rate, n_samps);  % win_limit_inner is the length of the *inner* time lag window, where data shoudl be excluded (always in one temporal direction).
if opt.isDirData
    pos = castFun(x);
elseif opt.longestDist
    pos = castFun(  round(   cumsum(  [0; sqrt((diff(x).^2) + (diff(y).^2))]  )   )    );
else
    pos = [castFun(x), castFun(y)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is Dan Manson's code for getting the exact nBSpikes in each window,   %
% hence pre-allocating the cross-corr arrays precisely.                      %
% TW: I have modifed this so that there can also be an 'inner window', for which data
% is excluded. This is how opt.offsetTimeWindow is implemented.
% 1a. For each A spikeInd, calculate first and last B spikes in the window 
%    and thus the number of B spikes.
[~,startOfWindowBSpkInd] = histc(spk_sampA - win_limit, spk_sampB);  % Don't be confused by the name, this is not a list of the first B spikes in each window, rather it is all the B spikes converted to the index: which spike A window (start) do I immediately follow?
startOfWindowBSpkInd = startOfWindowBSpkInd + 1;
[~,endOfWindowBSpkInd] = histc(spk_sampA + win_limit,spk_sampB);  % Likewise,this is the B spike train converted to index: 'which spike A window (end) do  I immediately precede?
startOfWindowBSpkInd(~startOfWindowBSpkInd) = 1;    %   I think that this line is unneccesary - after you +1, two lines above, there can not be 1s in the vector.
endOfWindowBSpkInd(~endOfWindowBSpkInd) = n_spksB;
nBSpikesInWin = endOfWindowBSpkInd - startOfWindowBSpkInd + 1;   % If you ever get confused about this code, run an AC, and look at [spk_sampA, startOfWindowBSpkInd, endOfWindowBSpkInd] all together - you can clearly see that for each spk_sampA, the other two indices tell you which are the first and the last spikes in that window. you 
% Inserted by TW: duplicate the above (1a), but calculate the inner window, spikes to be excluded.
% I also had to change the code a bit: I was getting zeros at *both* ends of the vectors xxxOfInnerWindowBSpkInd (I can't understand why this doesn't happen in the original code).
% Zeros at the beginning are correct (for InnerWindow), only those at the end need to be set to n_spksB.
% To fix, I added length(x) to the end of the histogram bins, then set values in xxxOfInnerWindowBSpkInd that were n_spksB+1, to n_spksB.
[~,startOfInnerWindowBSpkInd] = histc(spk_sampA - win_limit_inner, [spk_sampB; length(x)] );
startOfInnerWindowBSpkInd( startOfInnerWindowBSpkInd==n_spksB+1 ) = n_spksB;
startOfInnerWindowBSpkInd = startOfInnerWindowBSpkInd + 1;
[~,endOfInnerWindowBSpkInd] = histc(spk_sampA + win_limit_inner, [spk_sampB; (length(x)+win_limit_inner)] );
endOfInnerWindowBSpkInd( endOfInnerWindowBSpkInd==n_spksB+1 ) = n_spksB;
% startOfInnerWindowBSpkInd(~startOfInnerWindowBSpkInd) = 1;
% endOfInnerWindowBSpkInd(~endOfInnerWindowBSpkInd) = n_spksB;
nBSpikesInInnerWin = endOfInnerWindowBSpkInd - startOfInnerWindowBSpkInd + 1;
% Subtract the N spikes in the inner window from those across the whole window %
nBSpikesInWin = nBSpikesInWin - nBSpikesInInnerWin;
%1b. Work out offset inidices to be used when storing spike data
off_spike = cumsum([0 ; nBSpikesInWin]);

%1c. For each spike A, "clip" the start/end window index to stay within the 1:n_samps range of pos. 
startOfWindowPosInd = max(spk_sampA -win_limit, 1);
endOfWindowPosInd = min(spk_sampA + win_limit, n_samps);
nPosInWindow =  endOfWindowPosInd - startOfWindowPosInd + 1;
% Inserted by TW: duplicate the above (1c), but calculate the inner window, positions to be excluded.
startOfInnerWindowPosInd = max(spk_sampA - win_limit_inner, 1);
endOfInnerWindowPosInd = min(spk_sampA + win_limit_inner, n_samps)   +   1;
nPosInInnerWindow =  endOfInnerWindowPosInd - startOfInnerWindowPosInd - 1;
% Subtract N pos inner window from outer window, get the final numbers of spikes per window %
nPosInWindow = nPosInWindow - nPosInInnerWindow;
off_dwell = double(cumsum([0 ; nPosInWindow])); 

%1d. Pre-allocate dwell and spike arrays. %
p_disp = zeros(off_dwell(end), nDataDims+1, typeString);   % 'nDataDims + 1' is to leave a 3rd column, for offset time.
s_disp = zeros(off_spike(end), nDataDims+1, typeString);   %
% Finish DM Cde                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-allocate arrays for shuffled data %
if opt.posShuffleActive
    % Make the shuffled pos array, which is used to construct the shuffled spike maps. This is a 2D array, where each column is the entire disp dataset, but wrapped around a 
    % fixed time point, changing column to column. When you look at each row, therefore, you get the displacement when the data has been offset by offset:offset:trialLength values.
    % Then, below in the "1:n_spkA" loop, this will be used to assign 'shuffled' directions to spikes in the windowed correlogram.
    shufOffsets = int16( round( linspace( opt.posShuffleMinOffset*pos_sample_rate, n_samps-(opt.posShuffleMinOffset*pos_sample_rate), opt.nPosShuffles) ));
    shufInd = bsxfun(@plus, shufOffsets, int16((1:n_samps)') );    
    shufInd(shufInd > n_samps) = shufInd(shufInd > n_samps) - n_samps;
    shufInd = permute(shufInd, [1 3 2]);    % Switch nRepShuf to 3rd dim.
    if opt.isDirData || opt.longestDist
        shuffledPos = x(shufInd);
    else
    	shuffledPos = cat(2, x(shufInd), y(shufInd));
    end
    % Pre-allocate shuffled data arrays.
    p_disp_shuf = zeros(size(p_disp,1), nColInShufArray, opt.nPosShuffles, typeString);
    s_disp_shuf = zeros(size(s_disp,1), nColInShufArray, opt.nPosShuffles, typeString); 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct windowed displacement data %
for ii = 1:n_spksA;
    ti = spk_sampA(ii);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Pos samples %%%
    window_pos = ( [ (startOfWindowPosInd(ii):startOfInnerWindowPosInd(ii)), endOfInnerWindowPosInd(ii):endOfWindowPosInd(ii) ] )';    % 'window_pos' is an index into the vector of pos samples, from the time of the current spkA at 'ti' to the end of the time lag (or the end of the trial). How to fix (for spikes as well)? 
    p_disp((1:nPosInWindow(ii)) + off_dwell(ii), :) = [bsxfun(@minus, pos(window_pos,:), pos(ti,:)), castFun(window_pos-ti)]; 
    % pos offsets for shuffled positions %
    if opt.posShuffleActive
        pDispShufTemp =  bsxfun(@minus, shuffledPos(window_pos,:,:), shuffledPos(ti,:,:) );     
        pDispShufTime = bsxfun(@minus, shufInd(window_pos,1,:), repmat(shufInd(ti,1,:), [1, size(shuffledPos,2), 1]) );  % Mark when p_disp_shuf goes over the end-of-the trial wrapping point in the shuffle.
        outerWinLength = endOfWindowPosInd(ii) - startOfWindowPosInd(ii) + 1;                                            % (to do this, check if the time offset is longer than the window length).
        pDispShufTemp( abs(pDispShufTime)>outerWinLength ) = 32768;                                                      % Use 32768 as a marker, as can't use NaN - assume that real data will never have such a large offset.
        p_disp_shuf((1:nPosInWindow(ii)) + off_dwell(ii), 1:nDataDims, :) = pDispShufTemp;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Spikes %%%
    window_spk_ind = startOfWindowBSpkInd(ii):endOfWindowBSpkInd(ii);                          % Get the spikes (B) in the time window (units=pos samp N)
    if nBSpikesInInnerWin(ii)
        window_spk_ind = window_spk_ind(   ~ismember(window_spk_ind, startOfInnerWindowBSpkInd(ii):endOfInnerWindowBSpkInd(ii))   );
    end    
    window_spk = spk_sampB(  window_spk_ind  );
    s_disp((1:nBSpikesInWin(ii)) + off_spike(ii), :) = [bsxfun(@minus, pos(window_spk,:), pos(ti,:)), castFun(window_spk-ti)];     % To get the spike positions, index into the existing pos for the window_spk
    % Spike offsets for shuffled positions %
    if opt.posShuffleActive
        sDispShufTemp = bsxfun(@minus, shuffledPos(window_spk,:,:), shuffledPos(ti,:,:));  
        sDispShufTime = bsxfun(@minus, shufInd(window_spk,1,:), repmat(shufInd(ti,1,:), [1, size(shuffledPos,2), 1])); % Mark when p_disp_shuf goes over the end-of-the trial wrapping point in the shuffle.        
        sDispShufTemp( abs(sDispShufTime)>outerWinLength ) = 32768;                                                    % (to do this, check if the time offset is longer than the window length).
        s_disp_shuf((1:nBSpikesInWin(ii)) + off_spike(ii), 1:nDataDims, :) = sDispShufTemp;                            % Use 32768 as a marker, as can't use NaN - assume that real data will never have such a large offset.
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% If requested, collapse 2D xy to 1D displacement distance %%%
if opt.showAs1D
    p_disp(:,1) = int16( sqrt( single(p_disp(:,1)).^2 + single(p_disp(:,2)).^2 ) );
    s_disp(:,1) = int16( sqrt( single(s_disp(:,1)).^2 + single(s_disp(:,2)).^2 ) );
    p_disp(:,2) = int16(1);  % Keep the second column (formerly y) as a column of ones - 
    s_disp(:,2) = int16(1);  % this keeps everything consistent for the map making code.
    if opt.posShuffleActive
        p_disp_shuf(:,:,1) = int16( sqrt( single(p_disp_shuf(:,:,1)).^2 + single(p_disp_shuf(:,:,2)).^2 ) );   % These lines are very inefficient for long windows/large spike counts (memory expansion for 
        s_disp_shuf(:,:,1) = int16( sqrt( single(s_disp_shuf(:,:,1)).^2 + single(s_disp_shuf(:,:,2)).^2 ) );   % int16 -> single conversion is the problem, I assume). Could be fixed with a x,r->R LUT?
        p_disp_shuf(:,:,2) = int16(1); % Keep dimensions consistent
        s_disp_shuf(:,:,2) = int16(1); % with 2D
    end
end
% Another sort of related thing: for 1D longest distance (where 'position' alaways increases monotonically forwards in time), 
% collapse backwards and forwards in time displacements, by taking the absolute displacement value.
if opt.longestDist
    p_disp(:,1) = abs(p_disp(:,1));
    s_disp(:,1) = abs(s_disp(:,1));
    if opt.posShuffleActive
        p_disp_shuf = abs(p_disp_shuf);
        s_disp_shuf = abs(s_disp_shuf);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make maps. (1) basic set up %
% Bin maps %
p_disp(:,1:nMapDims) = p_disp(:,1:nMapDims) ./ binsize; % ROUND is automatic as p_disp is int16  %% TODO - this could be done with a LUT, and it would avoid the half-bin at zero.  (LUT is (1:nOffsets), LUT(i) = bin for offset i).
s_disp(:,1:nMapDims) = s_disp(:,1:nMapDims) ./ binsize; % ROUND is automatic as s_disp is int16  %% This is important, when using shuffled control, the half-bin at zero makes the first value of the AC smaller than the second. (undersampling leads to higher SD, hence lower Z)
% For ACCUMARRAY, need to convert x,y to r,c indices.
maxDispBin = max(max(abs((p_disp(:,1:nMapDims)))));
if opt.showAs1D  || opt.longestDist;   offsetLeft = 0;   else   offsetLeft = maxDispBin;   end
p_disp(:,1:nMapDims) = p_disp(:,1:nMapDims) + offsetLeft + 1;
s_disp(:,1:nMapDims) = s_disp(:,1:nMapDims) + offsetLeft + 1;
% Assign the X bin position vector - this is (1) required to interpret output, (2) acts as a indicator for AC output size %
binVect = -(offsetLeft*binsize)  :  binsize  :   (maxDispBin*binsize);
% Pre-assign output %
if nMapDims==1;   nMapCol=1;   else   nMapCol=length(binVect);  end        % If output rate map is 1D (radial offset or direction), this will be in the row dim.
if opt.offsetTimeWindow;   n3Dim=1;   else   n3Dim=size(T,1);   end        % Multiple time lags will be stacked in the 3rd dimension - multiple time lags are NOT compatible with an offset time-lag.
[posMap,spkMap,shuf.binMean,shuf.binSD,shuf.maps] = deal( nan( length(binVect), nMapCol, n3Dim ) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make maps (2). Create map for each specified Time Lag %
timeLagLims = T .* pos_sample_rate;
for ii=1:size(timeLagLims,1)
    % Filter displacement data so as to get only the appropriate time lags
    posIndForTimeLag = p_disp(:,end)>=-timeLagLims(ii,2) & p_disp(:,end)<=timeLagLims(ii,2);   % The col where time lag is stored will be =2 for dir data, =3 for x,y data.
    spkIndForTimeLag = s_disp(:,end)>=-timeLagLims(ii,2) & s_disp(:,end)<=timeLagLims(ii,2);   % 
    
    % Make Rate Map %
    posMap(:,:,ii) = accumarray(p_disp( posIndForTimeLag, 1:nMapDims),1,[length(binVect), nMapCol]);
    spkMap(:,:,ii) = accumarray(s_disp( spkIndForTimeLag, 1:nMapDims),1,[length(binVect), nMapCol]);
    
    % Make shuffled rate maps from spike counts %
    % %% TODO %% I still haven't implemented multiple time lag maps here. Old method wasn't going to work with 
    %            ACCUMARRAY - would have required creation of a copy (of a sub-set) of x_disp_shuf for each time
    %            lag, with accompanying memory problems. Best prob fix is create one single copy to act as a 
    %            backup (in a IF clause, only when multiple time lags in use), and then use xxxIndForTimeLag
    %            to set offset values in x_disp_shuf outside the time lag to 32768.
    if opt.posShuffleActive
        % Round and bin the data, using the same binning as for the real data %
        p_disp_shuf = (p_disp_shuf ./ binsize) + offsetLeft + 1; 
        s_disp_shuf = (s_disp_shuf ./ binsize) + offsetLeft + 1; 
        % Data outside the bounds of the real data needs to the set to a dummy bin value (which will later be removed from the map) %
        p_disp_shuf( p_disp_shuf<1 | p_disp_shuf>length(binVect) ) = length(binVect) + 1; 
        s_disp_shuf( s_disp_shuf<1 | s_disp_shuf>length(binVect) ) = length(binVect) + 1;
        % Add an explicit index for the nShufRep dimension - this goes in the last column of x_disp_shuf, which could be 2 or 3 depending on running mode %
        p_disp_shuf(:,end,:) = bsxfun(@plus,      zeros(  [size(p_disp_shuf,1), 1, opt.nPosShuffles]  ,typeString),       permute(castFun(1:opt.nPosShuffles),[1 3 2])     );
        s_disp_shuf(:,end,:) = bsxfun(@plus,      zeros(  [size(s_disp_shuf,1), 1, opt.nPosShuffles]  ,typeString),       permute(castFun(1:opt.nPosShuffles),[1 3 2])     );
        % Reshape the offset array to column vectors (one col per dim) %
        p_disp_shuf = reshape(   permute(p_disp_shuf, [1 3 2]),   [], nColInShufArray );
        s_disp_shuf = reshape(   permute(s_disp_shuf, [1 3 2]),   [], nColInShufArray );
        % Make the maps %
        if opt.showAs1D || opt.isDirData || opt.longestDist                         % The shuffled map dimensions will depend on the running mode -
            shufMapSize = [length(binVect)+1, opt.nPosShuffles];                    % this allows us to save an extra column in x_disp_shuf when it is 
        else                                                                        % not necesssary (1D maps - radial disp or Dir).
            shufMapSize = [length(binVect)+1, length(binVect)+1, opt.nPosShuffles]; %
        end                                                                         % 
        posMapShuf = accumarray(  p_disp_shuf, 1, shufMapSize);
        spkMapShuf = accumarray(  s_disp_shuf, 1, shufMapSize);
        if opt.showAs1D || opt.isDirData || opt.longestDist
            posMapShuf = permute(posMapShuf, [1 3 2]);   % For the cases where the map is only 1D, we now shift the
            spkMapShuf = permute(spkMapShuf, [1 3 2]);   % shuf reps to the 3rd dim, for consistency with 2D maps.
            posMapShuf = posMapShuf(1:end-1,:,:); % Remove the dummy bin for bad offsets
            spkMapShuf = spkMapShuf(1:end-1,:,:); %                  ''
        else
            posMapShuf = posMapShuf(1:end-1,1:end-1,:); % Remove the dummy bin for bad offsets
            spkMapShuf = spkMapShuf(1:end-1,1:end-1,:); % - this is for 2D maps
        end
        % Smooth the SHUFFLED RATE and SHUFFLED POS maps %
        if ~isempty(opt.boxcar)
            kern = ones(opt.boxcar,opt.boxcar^(nMapDims-1),1);  % The length of the kernel in the 2nd dim is 1 for 1D maps, opt.boxcar for 2D maps. Always length 1 in the 3rd dim, so we don't smooth over different shuffles.
            kern = kern ./ sum(kern(:));
            posMapShuf = imfilter(posMapShuf,kern,0,'conv');  % These maps will have edge condition inaccuracies (smoothing assumes zero outside image), 
            spkMapShuf = imfilter(spkMapShuf,kern,0,'conv');  % but this will be corrected for when making rate map.
        end
        % Make the SHUFFLED RATE maps %    
        rateMapShuf = (spkMapShuf ./ posMapShuf) .* pos_sample_rate;  % spk/pos, and convert to Hz.
        % Remove bins not visited in at least 50% shuffles %
        underVisInd = bsxfun(@times, sum(posMapShuf>0,3)<(opt.nPosShuffles/2), ones(1,1,opt.nPosShuffles));
        rateMapShuf( logical(underVisInd) ) = nan;
        % Get the various measures of dispersion from the shuffled maps %
        shuf.binMean(:,:,ii) = nanmean(rateMapShuf,3);
        shuf.binSD(:,:,ii) = nanstd(rateMapShuf,0,3);
        % Assign the shuffled maps to the output struct.
        shuf.maps = rateMapShuf;
    end
    
end

% Smooth real rate map %
if ~isempty(opt.boxcar)
    kern = ones(opt.boxcar,opt.boxcar^(nMapDims-1));  % The length of the kernel in the 2nd dim is 1 for 1D maps, opt.boxcar for 2D maps. Always length 1 in the 3rd dim, so we don't smooth over different shuffles.
    kern = kern ./ sum(kern(:));
    posMap = imfilter(posMap,kern,0,'conv');  % These maps will have edge condition inaccuracies (smoothing assumes zero outside image), 
    spkMap = imfilter(spkMap,kern,0,'conv');  % but this will be corrected for when making rate map.
end

% Filter for points of low dwell time %
posMap(posMap < opt.dwellThresh*pos_sample_rate) = nan;

% Make ratemap output %
rateMap = (spkMap./posMap) .* pos_sample_rate;   % Convert to Hz as we are mkaing the rate map.

% Normalise to shuffle: convert rate to S.D's of shuffle beyond shuffle mean (i.e. convert to Z-scores of shuffled population) %
if opt.posShuffleActive && opt.normaliseToShuf
    rateMap = (rateMap - shuf.binMean)  ./  shuf.binSD;
    rateMap( shuf.binSD==0 ) = nan;
end

rtn = rateMap;
binCoordsOut = double(binVect);

clear p_disp s_disp full_p_disp full_s_disp
if opt.posShuffleActive
    clear s_disp_shuf p_disp_shuf posShufInd shufInd shuffledDir
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outMap]=smoothFunc(inMap,opt)
% Smooth maps %
if opt.showAs1D || opt.isDirData || opt.longestDist
    % 1D data %
    outMap = smooth(inMap,opt.boxcar,'moving');
    
else
    % 2D Data %
    unVisInd = isnan( inMap );
    padMap = nan(  size(inMap)+opt.boxcar-1  );
    padMap(   (1:size(inMap,1)) + ((opt.boxcar-1)/2),   (1:size(inMap,2)) + ((opt.boxcar-1)/2)    ) = inMap;
    colMap = im2col(padMap,[opt.boxcar, opt.boxcar],'sliding');
    outMap = reshape(   nanmean(colMap,1),  size(inMap)  );
    outMap( unVisInd ) = nan;
    
end
















