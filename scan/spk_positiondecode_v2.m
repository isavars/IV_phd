function [rtn] = spk_bayesdecode(T, maps, posMap, spikeTimes, posSamp, binSize, varargin)
% Bayesian decode of position from spike train. (Assuming spikes follow a Poisson process with a mean described by the rate map).
% See Zhang, Ginzberg, Sejnowski, McNaughton, 
%
%       [rtn] = spk_positiondecode_v2(T, maps, posMap, spikeTimes, posSamp, binSize)
%
% Input 'binSize' is the number of pixels in each bin, assuming that you have passed the x and the y in pixels as units.
%
% rtn.errors = decodeErrors;
% rtn.predProbs = maxPVals;
% rtn.x = actualX;
% rtn.y = actualY;
% rtn.xPred = expX(mostLikelyBins);
% rtn.yPred = expY(mostLikelyBins);
% rtn.pVals = pVals;
% rtn.visInd = visInd;
% rtn.timeWinsWithoutSpikes = sum(timeWinsWithoutSpikes) / length(timeWinsWithoutSpikes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse Optional Input parameters %
prms.posSampRate = 50;
prms.usePriorPosition = 0;
prms.alphaForPriorPos = 1;
prms.alphaSpeedScaled = 0;
prms.nanPerWinLimit = 0.5;   % If the proportion of x or y in a time window T is greater than this, do not use (set x,y for whole time window to nan)
prms.speed = [];
prms.dirMode = 0;
prms.swrMode = 0;
if isstruct(varargin{1})
    optIn = varargin{1};
    f=fieldnames(optIn);
    for ii=1:length(f);   prms.(f{ii}) = optIn.(f{ii});   end
elseif ischar(varargin{1})
    for ii=1:2:length(varargin)
        prms.(varargin{ii}) = prms.(varargin{ii+1});
    end
end
% How we get the last time point, and number of bins, depends on if we are reconstructing real data or SWRs %
if ~prms.swrMode
    decodeEndTime = trialDur;
    nDecodeWin    = floor( decodeEndTime./T ); 
else
    nDecodeWin    = ceil(  max(cellfun(@max,spikeTimes)) ./ T   );
    decodeEndTime = nDecodeWin .* T;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Pre-filter spikes for when the position data is set to NaN %
% We use this filtering to make sure that any filtering applied to the position data (e.g. speed filtering, 
% edge re-scaling) is also applied to the spike train before reconstruction.                                                                             %
nanPosInd = find(any(isnan(posSamp,2)));
for ii=1:length(spikeTimes)
    spikeTimesInPosSamps=ceil(spikeTimes{ii} ./ posSampRate);
    spikeTimesInPosSamps(spikeTimesInPosSamps==0) = 1;
    spikesToRemove = ismember(spikeTimesInPosSamps,nanPosInd);
    spikeTimes{ii} = spikeTimes{ii}(~spikesToRemove);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Generate 2D matrix of actual observed firing during each decode time window. %
% Format is (1:nDecodeWin,1:nCell).
actualFiring = nan(nDecodeWin,length(maps));  % Preallocate
for ii=1:length(spikeTimes)
    temp=histc(spikeTimes{ii},0:T:decodeEndTime);
    actualFiring(:,ii) = temp(1:end-1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Generate 2D matrix of expected firing values in each spatial bin % 
% (i.e. rate map value * decode window duration. Format is (1:nSpatialBin,1:nCell) %
expFiring = nan(     numel(maps{1}),    length(maps)   );     % Preallocate - N Rows is number of bins in map.
% Loop through supplied maps, linearising each map into a column of the expected firing array %
for ii=1:length(maps)
    expFiring(:,ii) = maps{ii}(:);    
end
expFiring = expFiring .* T;               % Convert mean rate to expected spikes in decode time window duration.
posMap = posMap(:) ./ nansum(posMap(:));       % .. and we also need a linearised dwell time probability map.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4) Now calculate the decode probabilities %
expFiring3D    = repmat(  shiftdim(expFiring',-1)  ,[nDecodeWin 1 1]);    % First get the actual and expected spikes in a format such that they can be directly operated on together,
actualFiring3D = repmat(actualFiring,[1 1 size(expFiring,1)]);            % ..format is (actual,nCell,Expected), i.e. (1:nDecodeWin, 1:nCell, 1:nSpatialBin), for both.     
expFiring3D( isnan(expFiring3D) ) = 0;    %  NaNs come from unvisited spatial bins. Need to do this as FACTORIAL won't operate on NaNs. We will remove the probs for these spatial bins below.                            
pValsEachCell = (   (expFiring3D.^actualFiring3D)  ./  factorial(actualFiring3D)   ) .* exp( - expFiring3D );   % This is the key equation: for each cell, what is the probability of seeing observed spiking, assuming a poisson process with a mean determined by the rate map.
pVals = squeeze( prod( pValsEachCell, 2 ) );                    % Get the overall P(spk|pos) distributions for each time window, assuming independent poisson spiking for each cell, by multiplying the by-cell probabilities. See Zhang et al 1998, Eq 35. Format is (1:nDecodeWin, 1:nBin).
pVals = pVals .* repmat(posMap',size(actualFiring,1),1);        % Multiply by the overall position probability distribution (i.e. pos map). Probs for univisted bins become NaN here again. pVals then becomes P(spk|pos)*P(pos), which is equivalent to P(pos|spk), by Bayes rule. (If you ignore P(spk), see Zhang et al p)
pVals = pVals ./ repmat(  nansum(pVals,2), 1, size(pVals,2)  ); % pVals are scaled such that, in each time window, the p for all spatial bins sums to 1.
[maxPVals,mostLikelyBins] = nanmax(pVals,[],2);    % This tells you which is the most likely bin, and how likely you were to be there.
% Assign Output %
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
% (6) Now get the actual mean position in each decode window, and find the error (actual - decoded) %
%    (only if comparing decode to actual position, not necessary if decoding SWRs).
if ~prms.swrMode
    
    %%%%%% (6a) Get the mean actual position in each decode window. %%%%%%
    % When taking mean positions, need to account for a) posSampRate might be 46.875Hz, b) could be directional (circular) data.
    posSampTimes = linspace(1/prms.posSampRate, size(posSamp,1)/prms.posSampRate, size(posSamp,1));
    posSampAsDecodeWins = ceil(  posSampTimes./T  );
    if ~isDirMode
        actualX = accumarray( posSampAsDecodeWins,  posSamp(:,1),  [nDecodeWin, 1],  @nanmean);
        actualY = accumarray( posSampAsDecodeWins,  posSamp(:,2),  [nDecodeWin, 1],  @nanmean);
    elseif isDirMode
        actualY = accumarray( posSampAsDecodeWins,  circ_ang2rad(posSamp(:,2)),  [nDecodeWin, 1],  @circ_mean);
        actualY = circ_rad2ang( actualY );
        actualY(actualY<0) = actualY(actualY<0) + 360;
        actualX = ones(size(actualY));   % actualX is a dummy for directional data.
    end
    % Assess the proportions the raw pos data set to NaN in each decode window. Windows with excessive NaNs will not be used. %
    nanPerDecodeWin = accumarray( posSampAsDecodeWins,  isnan(posSamp,2),  [nDecodeWin, 1])  ./  accumarray( posSampAsDecodeWins,  ones(size(posSamp,1),1),  [nDecodeWin, 1]);  % The complex denominator is needed for when SR=46.875Hz
    tooManyNanPerWinInd = nanPerDecodeWin>=prms.nanPerWinLimit;
    
    %%%%%% (6b) Get the distance error between mean actual position and decoded position for each decode window %%%%%%
    % Make an x and y co-ord lookup, to be able to retrive error in 2D %
    [expX, expY] = meshgrid(1:size(maps{1},2),1:size(maps{1},1));
    expX = (expX.*binSize) - (binSize/2);          % Convert the x,y from bin counts to actual position
    expY = (expY.*binSize) - (binSize/2);          %  ..
    % Calculate the error.
    decodeErrors = sqrt(  (actualX' - expX(mostLikelyBins)).^2   +   (actualY' - expY(mostLikelyBins)).^2   );    % This still has the format (1:nWin)
    if isDirMode
        % For dir data, errors so far will be in the range 0-360. Errors > 180 need to be mapped onto range 0-180. 
        decodeErrors( decodeErrors>180 )  = abs(  decodeErrors( decodeErrors>180 ) - 360  );   % -360 gets the circular equivalent, this is in range -180-0, but we only want the absolute error.
    end
    maxPVals( tooManyNanPerWinInd )     = nan;  % Don't use errors from windows with too many NaNs in raw pos data.
    decodeErrors( tooManyNanPerWinInd ) = nan;
    
    %%%%%% (6c) Additionally, we can recover p(pos|spk) at each actual visited position %%%%%%%
    % (i.e. never mind the max pVal, at the what was the pVal at the *actual* position of the rat, in each decode window) %
    binnedY = ceil(actualY./binSize);                          % Convert the mean actual postion per decode window from pixels to bins.
    binnedX = ceil(actualX./binSize);                          %  ..
    binnedY(binnedY>size(maps{1},2)) = size(maps{1},2);        % These two lines here are to prevent rounding errors arising from CEIL occuring in the call to sub2ind, below.
    binnedX(binnedX>size(maps{1},1)) = size(maps{1},1);        %  ..
    actualPosLinInd = sub2ind(       size(maps{1}),    binnedY,     binnedX        );              % Then, convert actualX and Y (format (1,nWin)) to a linear index into the rate map ..
    actualPosLinInd = actualPosLinInd   +   0 : numel(posMap) : ((numel(posMap)*(nDecodeWin-1));   % and convert to one linear index into pVals (1:nDecodeWin, 1:nSpatialBin).
    pValsActualPos  = pVals( actualPosLinInd );
    
    %%%%% (6d) Assign output from this section %%%%
    rtn.errors = decodeErrors;
    rtn.x = actualX';   % Flip these round so that everything is a column vector
    rtn.y = actualY';   %
    rtn.xPred = expX(mostLikelyBins);
    rtn.yPred = expY(mostLikelyBins);
    rtn.pValsActualPos = pValsActualPos;
    rtn.perTimeWinInUse = 1 - nanPerDecodeWin';
    rtn.timeWinsInUse = sum(~tooManyNanPerWinInd) / nDecodeWin * 100;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (7) If requested, calculate 'coherence' of actual firing versus firing predicited from decode        %
% See Jackson & Redish, 2003, Network Comp Neural Sys. for ref of method (eq numbers from this paper)  %
% Also note that this function only returns 'incoherence', the difference between actual and predicted %
% firing. To go from there to 'coherence' is done use a population-level bootstrap.                    %
if prms.getCoherence
    %%%% (7a) First, get the 'activty packet' of the actual, observed firing (eq 9) %%%%%
    normTuningCurves = bsxfun(@rdivide, expFiring, nanmax(expFiring,[],1));             % 'normTuningCurves' are the (linearised) rate maps, scaled so max=1. Format is (1:nSpatialBin,1:nCell).
    normTuningCurves = permute(normTuningCurves,[3, 2, 1]);                             % Get the dimensions so as to multiply directly with actualFiring (i.e. 1, nCell, nSpatialBin).
    actPackActual = nansum(   bsxfun(@times, normTuningCurves, actualFiring),  2   ) ;  % Activity packet (actual): mulitply tuning curves by actual firing, i.e. nspikes, for each cell and decode window. Then sum across cells.
    actPackActual = bsxfun(@rdivide, actPackActual, nansum(normTuningCurves, 2) );      % Actvity packet is normalised by the sum of all tuning curves (across cells) (also eq 9).
    
    %%%%% (7b) Similarly, the 'expected' activty packet. %%%%% 
    % Instead of multiplying the tuning curves by the actual observed firing, here we multiply them by the firing that is predicted 
    % by the decoded position (by reading off the expected firing (variable 'expFiring') what this should be at each decoded position). (Eq 12)
    expFiringFromDecode = expFiring( mostLikelyBins, : );   % Get expected n spikes for each cell and each decode window, at decoded position. Read it off rate map (actually map/T for rate->spikes). Format = (1:nDecodeWin,1:nCell).
    actPackExp          = nansum(   bsxfun(@times, normTuningCurves, expFiringFromDecode),  2   );  % Activity packet (expected).
    actPackExp          = bsxfun(@rdivide, actPackExp, nansum(normTuningCurves, 2) );               % Actvity packet is normalised by the sum of all tuning curves (across cells).
    
    %%%%% (7c) Compare actual and expected actvity packets %%%%%
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
    






                                            










