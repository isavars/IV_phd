function [rtn] = spk_bayesdecode(T, maps, posMap, spikeTimes, posSamp, binSize, varargin)
% Bayesian decode of position from spike train. (Assuming spikes follow a Poisson process with a mean described by the rate map).
% See Zhang, Ginzberg, Sejnowski, McNaughton, 1998 J. Neurophysiol.
%
%       [rtn] = spk_bayesdecode(T, maps, posMap, spikeTimes, posSamp, binSize)
%       [rtn] = spk_bayesdecode(T, maps, posMap, spikeTimes, posSamp, binSize, optionalInputStruct)
%       [rtn] = spk_bayesdecode(T, maps, posMap, spikeTimes, posSamp, binSize, 'inputName', inputVal, .. etc .. )
%
% Inputs: T          - Decode time window, in sec
%         maps       - cell array {1:nCell} of binned rate maps.
%         posMap     - binned position map, i.e. dwell time map.
%         spikeTimes - cell array, {1:nCell}, of spike trains as spike times, in sec.
%         posSamp    - raw pos data in DACQ pos samples, either XY (1:nPosSamp,2), or dir (1:nPosSamp,2).
%         binSize    - size of spatial bin, in the same units as 'posSamp'.
% 
% Optional inputs/analysis parameters (supply as struct or " ,'fieldname',value, " comma-separated list) :
%
% prms.posSampRate = 50;
% prms.usePriorPosition = 0;
% prms.alphaForPriorPos = 1;
% prms.alphaSpeedScaled = 0;
% prms.nanPerWinLimit = 0.5;   % If the proportion of x or y in a time window T is greater than this, do not use (set x,y for whole time window to nan)
% prms.speed = [];
% prms.dirMode = 0;
% prms.swrMode = 0;
% prms.getCoherence = 0;
% prms.overlappingT = 0; %y/n
% prms.winShift     = 0;
%
% Output structure:
%
% rtn.pValsBest  - Max pVal for each decode window.
% rtn.posPredInd - This is the index into the linearised rate map that gives the predicited position - contrast with xPred and yPred, below, which give actual predicited X and Y.
% rtn.pVals      - This is the full (1:nDecodeWin, 1:nBin) array of decode probs.
% rtn.errors     - distance error (decoded-actual) for each window.
% rtn.x          - Actual X position (mean for decode window)
% rtn.y          - Actual Y (or Dir) position (mean for decode window)
% rtn.xPred      - Decoded X position.
% rtn.yPred      - Decoded Y (or Dir) position.
% rtn.pValsActualPos  - probability from decode of rat being at its actual position
% rtn.incoherence     - incoherence score (difference between actual and expected actvity packets) for each decode window.
% rtn.perTimeWinInUse - for each decode window, how much raw pos data was set to NaN (from pre-existing filtering procedures).
% rtn.timeWinsInUse   - how many decode windows were actaully used (not too man nans)?




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse Optional Input parameters %
prms.posSampRate = 50;
prms.usePriorPosition = 0;
prms.alphaForPriorPos = 1;
prms.alphaSpeedScaled = 0;
prms.nanPerWinLimit   = 0.5;   % If the proportion of x or y in a time window T is greater than this, do not use (set x,y for whole time window to nan)
prms.speed            = [];
prms.dirMode          = 0;
prms.swrMode          = 0;
prms.useFullTrialSWR  = 0;
prms.trialLength      = nan;
prms.getErrors        = 1;
prms.getCoherence     = 0;
prms.overlappingT     = 0;
prms.winShift         = 0;
prms.useGPU           = 0;
%LM bugfix
if nargin>6 && isstruct(varargin{1})
    optIn = varargin{1};
    f=fieldnames(optIn);
    for ii=1:length(f);   prms.(f{ii}) = optIn.(f{ii});   end
elseif nargin>6 && ischar(varargin{1})
    for ii=1:2:length(varargin)
        prms.(varargin{ii}) = varargin{ii+1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-processing of the data, calculation of decode windows etc.
% Check if T is divisible by the shifting window parameter (if supplied). Bin construction assumes it is.
if prms.overlappingT 
    if ~(mod(T,prms.winShift)==0) 
        prms.winShift = T/2;
        warning('Decoding window length not divisible by window shift. Set to %.3f seconds (T/2) instead',T/2);
    end
    winShiftPerT  = T / prms.winShift;
end
% Get the start point of the first bin, and number of bins. How this works depends on if we are reconstructing real 
% data or just prediciting position from SWRs.
if ~prms.swrMode || prms.useFullTrialSWR
    % Decode the whole trial %
    decodeStartTime = 0;
    if ~isempty(posSamp)
        nDecodeWin = floor( (size(posSamp,1)/prms.posSampRate) ./ T );
    elseif ~isnan( prms.trialLength )
        nDecodeWin = floor( prms.trialLength ./ T );
    else
        error('I''m stuck: I don''t know how long the trial is');
    end
    % If overlappingT, the actual number of bins will be almost T*winShift, but not quite, as the 
    % final overlaps shouldn't go past the end of the trial. 
    if prms.overlappingT   
        nDecodeWin = (round((T/prms.winShift)) * nDecodeWin) - (winShiftPerT-1);
    end
else
    % SWR Mode - just decode from the first to the last spike passed as input %
    maxSpike       = max(cellfun(@max,spikeTimes(~cellfun('isempty',spikeTimes))));
    minSpike       = min(cellfun(@min,spikeTimes(~cellfun('isempty',spikeTimes))));
    if prms.overlappingT;   
        nDecodeWin = ceil( (maxSpike - minSpike) / prms.winShift ) - (winShiftPerT-1);
    else
        nDecodeWin = ceil( (maxSpike - minSpike) / T );
    end
    if nDecodeWin == 0; nDecodeWin = 1; end % protect if only 1 spike in window (should probably escape earlier in that case)
    decodeStartTime = minSpike;
end
posSamp = double( posSamp );   % If calling from scan, pos data is integers, make double.
% Force the cell arrays for spikes and maps into the correct orientation, otherwise the 
% 'cell2mat' operations below to get their contents into single arrays won't work.
if size(spikeTimes,2) > size(spikeTimes,1)
    spikeTimes = spikeTimes';
end
if size(maps,2) < size(maps,1)
    maps = maps';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Pre-filter spikes for when the position data is set to NaN %
% We use this filtering to make sure that any filtering applied to the position data (e.g. speed filtering, 
% edge re-scaling) is also applied to the spike train before reconstruction.
% This is not necessary for SWR mode: here we assume the spikes are 'offline' firing, and therefore any filtering
% applied to 'online' positions is not relevant. 
if ~prms.swrMode   
    nanPosInd = find(  any(   isnan(posSamp),  2   )  );
    for ii=1:length(spikeTimes)
        spikeTimesInPosSamps=ceil(spikeTimes{ii} .* prms.posSampRate);  % Bug fix here, TW, 2017-07-17: .*prms.posSampRate was ./prms.posSampRate, meaning that code was potentially broken for nans in pos data (all spikes removed if nans early in trial).
        spikeTimesInPosSamps(spikeTimesInPosSamps==0) = 1;
        spikesToRemove = ismember(spikeTimesInPosSamps,nanPosInd);
        spikeTimes{ii} = spikeTimes{ii}(~spikesToRemove);
        % bug fix for cells with no spikes - there is a difference in format 
        % between cells that fire no spikes in whole trial or just not in 
        % decoding window
        if size(spikeTimes{ii},2) == 0
            spikeTimes{ii} = zeros(0,1); 
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Generate 2D matrix of actual observed firing during each decode time window. %
% Format is (1:nDecodeWin,1:nCell).
%%% TW edit: new, faster method of creating 'actualfiring' spk histogram %
% Bin the spikes into time windows: bins are either T, for 'normal' operation, or prms.winShift, for sliding bins %
if prms.overlappingT;   t = prms.winShift;   else    t = T;   end
spkVect             = cell2mat( spikeTimes );
spkVect             = ceil( (spkVect-decodeStartTime) ./ t );
spkVect(spkVect==0) = 1;                                         % bin values will be used as an index, so 0 isn't allowed.
% Make another index, of which spike belongs to which cell, (to build the (nWin, nCell) array of spk counts). %
indCell = cellfun( @(x) ones(length(x),1), spikeTimes, 'UniformOutput', 0 );
for ii=1:length(indCell);   indCell{ii}=indCell{ii}.*ii;   end
indVect             = cell2mat(indCell);
% Now build the spk array: if overlapping windows, first build the array of prms.winShift bins, 
% then convolve with a kernel of winShiftPerT.
if ~prms.overlappingT
    actualFiring    = accumarray( [spkVect, indVect], 1, [nDecodeWin, length(spikeTimes)]  );
else
    spkArr          = accumarray( [spkVect, indVect], 1, [nDecodeWin+(winShiftPerT-1), length(spikeTimes)]  );
    k               = ones( ceil(winShiftPerT), 1 ); %% num format error when e.g. 0.3/0.1
    actualFiring    = conv2( spkArr, k, 'valid' );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Generate 2D matrix of expected firing values in each spatial bin % 
% (i.e. rate map value * decode window duration. Format is (1:nSpatialBin,1:nCell) %
if size(maps,1)>size(maps,2);   maps=maps';   end             % Force rate map cell array into column cell-vector ..
mapsLin   = cellfun( @(x) x(:), maps, 'UniformOutput', 0 );   % .. and rate maps into row vectors, so that the 'cell2mat' below always gives a 2D, (spatBin,cell), array. 
expFiring = cell2mat( mapsLin );
expFiring = expFiring .* T;                       % Convert mean rate to expected spikes in decode time window duration.
expFiring( isnan(expFiring) ) = 0;                %  NaNs come from unvisited spatial bins. Need to do this as FACTORIAL won't operate on NaNs. We will remove the probs for these spatial bins below.         
posMap    = posMap(:) ./ nansum(posMap(:));       % .. and we also need a linearised dwell time probability map.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4) Now calculate the decode probabilities .. %
% The trick here is to vectorise the whole calculation, on a (nTimeBin, nCell, nSpatBin) array. The key equation
% is coded as a function handle, and applied to all values using bsxfun.
% Can also use the GPU: principle is similar, use ARRAYFUN instead of BSXFUN but as ARRAYFUN on a GPU doesn't
% support FACTORIAL, we need to calulate this separately and feed it as a separate argument. 
expFiring     = shiftdim(expFiring', -1);  
function [pForCell] = eqFunGPU( actF_g, actFF_g, expF_g )
    pForCell = ( (expF_g.^actF_g) ./ actFF_g ) .* exp( - expF_g );
end
if ~prms.useGPU
    eqFun         = @(expF,actF) ( (expF.^actF) ./ factorial(actF) ) .* exp( - expF );  % Equation defining P(spk|pos), for each cell assuming poisson spiking. See Zhang et al 1998, Eq 35.
    pValsEachCell = bsxfun( eqFun, expFiring, actualFiring );
    pVals         = prod( pValsEachCell, 2 );                                           % Get the overall P(spk|pos) distributions, by multiplying the by-cell probabilities.
else
    actF_g        = gpuArray( actualFiring );
    expF_g        = gpuArray( expFiring );
    actFF_g       = factorial( actF_g );
    pValsEachCell_g = arrayfun( @eqFunGPU, actF_g, actFF_g, expF_g );
    pVals_g         = prod( pValsEachCell_g, 2 );
    pVals           = gather( pVals_g );
    
end
pVals = squeeze( pVals );
pVals = bsxfun( @times, pVals, posMap' );                % Multiply by the overall position probability distribution (i.e. pos map). Probs for univisted bins become NaN here again. pVals then becomes P(spk|pos)*P(pos), which is equivalent to P(pos|spk), by Bayes rule. (If you ignore P(spk), see Zhang et al p)
pVals = bsxfun( @rdivide, pVals, nansum(pVals,2));       % pVals are scaled such that, in each time window, the p for all spatial bins sums to 1.
[maxPVals,mostLikelyBins] = nanmax(pVals,[],2);          % This tells you which is the most likely bin, and how likely you were to be there.

rtn.pValsBest  = maxPVals;
rtn.posPredInd = mostLikelyBins; % This is the index into the linearised rate map that gives the predicited position - contrast with xPred and yPred, below, which give actual predicited X and Y.
rtn.pVals      = pVals;          % This is the full (1:nDecodeWin, 1:nBin) array of decode probs.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5) Combine 'online' position estimate with prior probability based on previous position %
%%% TODO - this still needs some debugging, have left this for now as this is rarely used %
% if prms.usePriorPosition
%     % If required, get the mean speed per bin:
%     speedByWin = nanmean(speedTemp);
%     % Actual speeds need to be linearly re-scaled into a range of alpha values. See Zhang et al 1998.
%     if scaledAlpha;   
%         alphaVals = (speedByWin .* alphaSpeedScale(2)) + alphaSpeedScale(1);
%     end  
%     
%     % Now need to cycle through time windows one-by-one, getting estimate on the basis of current ML and prior %
%     pTempCurr = zeros(mapSize);
%     pTempPrior = zeros(mapSize);
%     pValsPrior = zeros(size(pVals));    % To store all prior P dists, for debugging purposes rather than for calculation.
%     [X,Y] = meshgrid(1:mapSize(2),1:mapSize(1));   % Grid for construction of 2-D gaussian .. 
%     X = X.*binSize;      Y = Y.*binSize;          %  .. convert the x,y from bin counts to actual position (in pix)
%     X = X-(binSize/2);   Y = Y-(binSize/2);       %  ..
%     for ii=1:nDecodeWin
%         pTempCurr(:)=0;  pTempPrior(:)=0;             % Reset prob dist variables to all zero.
%         pTempCurr = pVals(ii,:);                      % Get the ML (spike-based) p distribution for the current time window. At this point, format changes from 1:nBin vector to 2D spatial map. 
%         if ii==1
%             pTempPrior=ones(size(pTempPrior));     % In first time bin, no info on prior, so flat distribution.
%         else
%             prevX=expX(predPosTemp);   prevY=expY(predPosTemp);   % Get the previous estimated position.
%             if scaledAlpha;   
%                 alpha=alphaVals(ii-1);   
%             else
%                 alpha=alphaVals;
%             end
%             exponent = (    ((X-prevX).^2)  +  ((Y-prevY).^2)    )    ./   (2*alpha*alpha); % 2D gaussian centred on previous estimated postion, (optional: with
%             pTempPrior = (exp(-exponent));                                                  % standard deviation scaled by speed in previous time window).               
%         end
%         pTempPrior = pTempPrior ./ nansum(pTempPrior(:));   % scaled to become an actual prob distribution that sums to 1 over the whole box.
%         pCurrTimesPrior = pTempCurr.*pTempPrior;
%         pCurrTimesPrior = pCurrTimesPrior ./ nansum(pCurrTimesPrior(:));
%         pVals(ii,:) = pCurrTimesPrior(visInd);             % New predicted position, incorporating prior prediction, is stored in variable 'pVals
%         pValsPrior(ii,:) = pTempPrior(visInd);   % Store all prior P dists, for debugging purposes rather than for calculation.
%         [~,predPosTemp] = max(pVals(ii,:));      % This is the current predicted position, will be needed to set the prior in next time window.
%     end
%           rtn.speed = speedByWin';   % Flip these round so that everything is a column vector
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6) .. and work out what are the positions associated with the maximum probabilities                %
% Make an x and y co-ord lookup, to be able to retrive error in 2D %
[expX, expY] = meshgrid(1:size(maps{1},2),1:size(maps{1},1));
expX = (expX.*binSize) - (binSize/2);          % Convert the x,y from bin counts to actual position
expY = (expY.*binSize) - (binSize/2);          %  ..
rtn.xPred = expX(mostLikelyBins);
rtn.yPred = expY(mostLikelyBins);       % If called with prms.dirMode=1, 'expY' (actaully dir) is still in degrees, as the user gives 'binSize' in degrees. xPred will be a dummy.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (7) Now get the actual mean position in each decode window, and find the error (actual - decoded) %
%    (only if comparing decode to actual position, not necessary if decoding SWRs).
if prms.getErrors && ~prms.swrMode  % Do double check - setting 'swrMode' to 1 should automatically mean no error calculation.
    
    %%%%%% (6a) Get the mean actual position in each decode window. %%%%%%
    % When taking mean positions, need to account for a) posSampRate might be 46.875Hz, b) could be directional (circular) data.
    posSampTimes = linspace(1/prms.posSampRate, size(posSamp,1)/prms.posSampRate, size(posSamp,1))';
    posSampAsDecodeWins = ceil(  posSampTimes./T  );
    if ~prms.dirMode
        actualX = accumarray( posSampAsDecodeWins,  posSamp(:,1),  [nDecodeWin, 1],  @nanmean);
        actualY = accumarray( posSampAsDecodeWins,  posSamp(:,2),  [nDecodeWin, 1],  @nanmean);
    else
        actualY = circ_rad2ang( accumarray( posSampAsDecodeWins,  circ_ang2rad(posSamp(:,1)),  [nDecodeWin, 1],  @circ_mean) );
        actualY(actualY<0) = actualY(actualY<0) + 360;  
        actualX = ones(size(actualY));   % actualX is a dummy for directional data.
    end
    % Assess the proportions the raw pos data set to NaN in each decode window. Windows with excessive NaNs will not be used. %
    nanPerDecodeWin = accumarray( posSampAsDecodeWins,  any(isnan(posSamp),2),  [nDecodeWin, 1])  ./  accumarray( posSampAsDecodeWins,  ones(size(posSamp,1),1),  [nDecodeWin, 1]);  % The complex denominator is needed for when SR=46.875Hz
    tooManyNanPerWinInd = nanPerDecodeWin>=prms.nanPerWinLimit;
    actualX( tooManyNanPerWinInd ) = nan;   actualY( tooManyNanPerWinInd ) = nan;   % Set all pos data for windows with too many NaNs to NaN.
    
    %%%%%% (6b) Get the distance error between mean actual position and decoded position for each decode window %%%%%%
    if ~prms.dirMode
        %%% LM COMMENT: Not sure why 'actualX' (and ...Y) were set to be transposed - this makes error = nWin:nWin which makes no sense, or does it? 
        decodeErrors = sqrt(  (actualX - expX(mostLikelyBins)).^2   +   (actualY - expY(mostLikelyBins)).^2   );    % For 2D data error is pythagoanen distance. This still has the format (1:nWin)
%         decodeErrors = sqrt(  (actualX' - expX(mostLikelyBins)).^2   +   (actualY' - expY(mostLikelyBins)).^2   );    % For 2D data error is pythagoanen distance. This still has the format (1:nWin)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        decodeErrors = circ_rad2ang(   circ_dist(circ_ang2rad(actualY), circ_ang2rad( expY(mostLikelyBins) ) )   );   % For direction data smallest distance across circle. The second arg needs converting to radians as the user will supply as degrees.
    end
    decodeErrors( tooManyNanPerWinInd ) = nan;     % Set all pos data for windows with too many NaNs to NaN.
    
    %%%%%% (6c) Additionally, we can recover p(pos|spk) at each actual visited position %%%%%%%
    % (i.e. never mind the max pVal, at the what was the pVal at the *actual* position of the rat, in each decode window) %
    binnedY = ceil(actualY./binSize);                   % Convert the mean actual postion per decode window from pixels to bins.
    binnedX = ceil(actualX./binSize);                   %  ..
    binnedY(binnedY>size(maps{1},1)) = size(maps{1},1);        % The 4 lines following here are to prevent rounding errors arising from CEIL occuring in the call to sub2ind, below.
    binnedX(binnedX>size(maps{1},2)) = size(maps{1},2);        %  .. also note the inversion of the r,c indices for the rate map, to match x,y ..
    binnedY(binnedY==0)              = 1;                      %  ..
    binnedX(binnedX==0)              = 1;                      %  ..
    actualPosLinInd = sub2ind(       size(maps{1}),    binnedY,     binnedX        );               % Then, convert actualX and Y (format (1,nWin)) to a linear index into the rate map ..
    actualPosLinInd = actualPosLinInd   +   (0 : numel(posMap) : (numel(posMap)*(nDecodeWin-1)))';   % and convert to one linear index into pVals (1:nDecodeWin, 1:nSpatialBin).
    posWinIsNan     = isnan( actualPosLinInd );  % Need to protect for the case where 
    actualPosLinInd(posWinIsNan) = 1;            % whole pos win is nan .. 
    pValsActualPos  = pVals( actualPosLinInd );
    pValsActualPos( posWinIsNan ) = NaN;         % .. set the actual PVals to NaN here.
    
    %%%%% (6d) Assign output from this section %%%%
    rtn.errors = decodeErrors;
    rtn.x = actualX;
    rtn.y = actualY;
    rtn.pValsActualPos = pValsActualPos;
    rtn.perTimeWinInUse = 1 - nanPerDecodeWin;
    rtn.timeWinsInUse = sum(~tooManyNanPerWinInd) / nDecodeWin * 100;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (8) If requested, calculate 'coherence' of actual firing versus firing predicited from decode        %
% See Jackson & Redish, 2003, Network Comp Neural Sys. for ref of method (eq numbers from this paper)  %
% Also note that this function only returns 'incoherence', the difference between actual and predicted %
% firing. To go from there to 'coherence' is done use a population-level bootstrap.                    %
if prms.getCoherence
    
    %%%% (7a) First, we need to make 'tuning curves' (i.e. (linearised) rate maps, scaled so max=1).
    normTuningCurves = bsxfun(@rdivide, expFiring, nanmax(expFiring,[],1));             % Format is (1:nSpatialBin,1:nCell).
    normTuningCurves = permute(normTuningCurves,[3, 2, 1]);                             % Get the dimensions so as to multiply directly with actualFiring (i.e. 1, nCell, nSpatialBin).
    
    %%%% (7b) Then, get the 'activty packet' of the actual, observed firing (eq 9) %%%%%
    actPackActual = nansum(   bsxfun(@times, normTuningCurves, actualFiring),  2   ) ;  % Activity packet (actual): mulitply tuning curves by actual firing, i.e. nspikes, for each cell and decode window. Then sum across cells.
    actPackActual = bsxfun(@rdivide, actPackActual, nansum(normTuningCurves, 2) );      % Actvity packet is normalised by the sum of all tuning curves (across cells) (also eq 9).
    
    %%%%% (7c) Similarly, the 'expected' activty packet. %%%%% 
    % Instead of multiplying the tuning curves by the actual observed firing, here we multiply them by the firing that is predicted 
    % by the decoded position (by reading off the expected firing (variable 'expFiring') what this should be at each decoded position). (Eq 12)
    expFiringFromDecode = expFiring( mostLikelyBins, : );   % Get expected n spikes for each cell and each decode window, at decoded position. Read it off rate map (actually map/T for rate->spikes). Format = (1:nDecodeWin,1:nCell).
    actPackExp          = nansum(   bsxfun(@times, normTuningCurves, expFiringFromDecode),  2   );  % Activity packet (expected).
    actPackExp          = bsxfun(@rdivide, actPackExp, nansum(normTuningCurves, 2) );               % Actvity packet is normalised by the sum of all tuning curves (across cells).
    
    %%%%% (7d) Compare actual and expected actvity packets %%%%%
    % Use the variance-of-difference method outlined in Johnson, Seeland & Redish, 2005, Hippocampus, rather than the RMSE method in Jackson et al 2003 % 
    % Note that both of these measure only the *difference* between activity packets %
    rtn.incoherence = nanvar( actPackActual-actPackExp, 0, 3 ) ./ nansum( actPackExp, 3 );   % Johnson et al 2005 Eq 6.    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------- Debug (plot probs, real and extimated positions) -------------------------------------------------%
if 0
    debugStart = 1;
    nRows=8;
    nCols=7;
    winInd = debugStart:debugStart+((nRows*nCols)-1);
    figure;   axCount=1;
    for ii=winInd
        subplot(nRows,nCols,axCount);  axCount=axCount+1;
        a = zeros(size(maps{1}));
        a(visInd) = pVals(ii,:);
        imagesc(a);
        hold on
        plot(actualX(ii)/10,actualY(ii)/10,'m*');
        plot(expX(mostLikelyBins(ii))/10, expY(mostLikelyBins(ii))/10, 'k+');
        axis off;   axis square;
        title(num2str(round(decodeErrors(ii))));
    end
    if usePriorPosition
        figure;   axCount=1;
        for ii=winInd
            subplot(nRows,nCols,axCount);  axCount=axCount+1;
            a = zeros(size(maps{1}));
            a(visInd) = pValsPrior(ii,:);
            imagesc(a,[0 5/(numel(maps{1}))]);    axis off;   axis square;
        end
    end
    figure; hist(decodeErrors,5:10:355);
    hold on
    p=prctile(decodeErrors,[50,75,90,95]);
    for jj=1:length(p);  plot([p(jj) p(jj)], get(gca,'ylim'), 'g-');   end
    m=mean(decodeErrors);
    plot([m m], get(gca,'ylim'), 'r-');
    set(gca,'xlim',[0 350]);
end
    

end




                                            










