function [rtn] = spk_positiondecode(decodeMethod,T,maps,posMap,spikeTimes,xSamp,ySamp,speed,trialDur,posSampRate,binSize,usePriorPosition,varargin)

% HELP
%       [rtn] = spk_positiondecode(decodeMethod,T,maps,posMap,spikeTimes,xSamp,ySamp,trialDur,posSampRate,binSize,usePriorPositionFlag)
%       [rtn] = spk_positiondecode(decodeMethod,T,maps,posMap,spikeTimes,xSamp,ySamp,trialDur,posSampRate,binSize,usePriorPositionFlag,alpha)
%       [rtn] = spk_positiondecode(decodeMethod,T,maps,posMap,spikeTimes,xSamp,ySamp,trialDur,posSampRate,binSize,usePriorPositionFlag,speed,alphaSpeedScaling)
%
% For dir data:
%
%       [rtn] = spk_positiondecode(decodeMethod,T,maps,posMap,spikeTimes,[],ySamp,trialDur,posSampRate,binSize,usePriorPositionFlag)
%       [rtn] = spk_positiondecode(decodeMethod,T,maps,posMap,spikeTimes,[],ySamp,trialDur,posSampRate,binSize,usePriorPositionFlag,smoothPValsFlag,kernLength)
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



%------------------------------------------------ Parse Input ------------------------------------------------------------------%
% Get the alpha for the prior 2D gaussian. Can either be fixed, or scale with speed, depening on user input %
scaledAlpha=0;
if usePriorPosition
    if length(varargin{1})==1
        scaledAlpha = 0;
        alphaVals = varargin{1};     % Fixed
    elseif length(varargin{1})==2
        scaledAlpha = 1;
        alphaSpeedScale = varargin{1};   % Scaling
    end
end

% Check for xSamp = []  -  this means that the input data is directional, we need to produce a dummy set of y data. %
isDirMode = 0;
smoothPMode = 0;
if isempty(xSamp)
    xSamp = ones(size(ySamp));
    if ~isempty(varargin) && varargin{1}
        smoothPMode=1;
    end
    isDirMode = 1;
end



%----------------------------- Pre-filter spikes for when the position data is set to NaN --------------------------------------%
% We can use this filtering to make sure that any filtering applied to the map (e.g. speed filtering, edge re-scaling) is also  %
% applied to the spike train before reconstruction.                                                                             %
xSamp( isnan(ySamp) ) = NaN;    % From here onwards, xSamp is used canonical source of NaN positions.
nanPosInd = find(isnan(xSamp));
for ii=1:length(spikeTimes)
    spikeTimesInPosSamps=ceil(spikeTimes{ii} ./ posSampRate);
    spikeTimesInPosSamps(spikeTimesInPosSamps==0) = 1;
    spikesToRemove = ismember(spikeTimesInPosSamps,nanPosInd);
    spikeTimes{ii} = spikeTimes{ii}(~spikesToRemove);
end



%-------------------------------- Mean the position and speed data into time windows -------------------------------------------%
% If sample rate=50, simple reshape. If DACQ1 (46.875Hz), a bit more long-winded %
nanPerTLimit = 0.5;   % If the proportion of x or y in a time window T is greater than this, do not use (set x,y for whole time window to nan)
nWin = floor(trialDur/T);
if posSampRate==50 && rem(T,1/posSampRate)==0
    % In the case when the sample rate is 50, and the time window is divisible by 1/50, simple reshape required, but may need to cut off 'overhanging' data that doesn't make a full time window % 
    xSamp=xSamp(1:round(nWin*T*posSampRate));  % Using 'round' as I once got some dodgy matlab round-off errors here!
    ySamp=ySamp(1:round(nWin*T*posSampRate));
    speed=speed(1:nWin*T*posSampRate);
    xTemp=reshape(xSamp,round(posSampRate*T),nWin);
    yTemp=reshape(ySamp,round(posSampRate*T),nWin);
    speedTemp = reshape(speed,posSampRate*T,nWin);   
else  
    % If sample rate=46.875, or the time window and sample rate are otherwise mis-matched, need to get mean pos data per T window using for looped index %
    posSampTimes = linspace(1/posSampRate, length(xSamp)/posSampRate, length(xSamp));
    posSampTimeWins = ceil(  posSampTimes./T  );
    xTemp=nan(ceil(posSampRate*T),nWin);
    yTemp=nan(ceil(posSampRate*T),nWin);
    if scaledAlpha;   speedTemp=nan(ceil(posSampRate*T),nWin);   end
    for ii=1:nWin
        ind=posSampTimeWins==ii;
        xTemp(1:sum(ind),ii) = xSamp(ind);
        yTemp(1:sum(ind),ii) = ySamp(ind);
        speedTemp(1:sum(ind),ii) = speed(ind);
    end
end
% Now that the X and Y data are reshaped into windows, we can take the means for each window.
if ~isDirMode
    actualX = nanmean(xTemp,  1  );
    actualY = nanmean(yTemp,  1  );
elseif isDirMode
    actualY = circ_rad2ang( circ_mean( circ_ang2rad(yTemp),  [],  1  ) );
    actualX = ones(size(actualY));
    actualY(actualY<0) = actualY(actualY<0) + 360;
end
speedByWin = nanmean(speedTemp);
% Windows which have too many nans in raw data are set to be NaN for the whole window.
propNanPerT = sum( isnan(xTemp), 1 ) ./ repmat(size(xTemp,1),1,size(xTemp,2));
actualX( propNanPerT >= nanPerTLimit ) = nan;
actualY( propNanPerT >= nanPerTLimit ) = nan;
speedByWin( propNanPerT >= nanPerTLimit ) = nan;
% Actual speeds need to be linearly re-scaled into a range of alpha values. See Zhang et al 1998.
if scaledAlpha;   
    alphaVals = (speedByWin .* alphaSpeedScale(2)) + alphaSpeedScale(1);
end  



%-------------------------------- Generate 2D matrix of actual observed firing during each time window ----------------------------------------------%
% Format is (1:nWin,1:nCell), where nWin = trial_dur/T.
actualFiring = nan(nWin,length(maps));  % Preallocate
for ii=1:length(spikeTimes)
    temp=histc(spikeTimes{ii},0:T:trialDur);
    actualFiring(:,ii) = temp(1:end-1);
end



%---------------------------- Generate 2D matrix of expected firing values, with format (1:nBin,1:nCell) --------------------------------------------%
visInd = ~isnan(maps{1});
expFiring = nan(     sum(visInd(:)),    length(maps)   );     % Preallocate - N Rows is number of visited bins, rather possible total in env.
% Loop through supplied maps, linearising each map into a column of the expected firing array %
for ii=1:length(maps)
    expFiring(:,ii) = maps{ii}(visInd);        % Lineraise using visInd.
    expFiring(:,ii) = expFiring(:,ii) .* T;    % Multiply by T to convert from Hz to expected spikes in T (time window)
end
% We also need an x and y co-ord lookup, to be able to retrive error in 2D %
[expX, expY] = meshgrid(1:size(maps{1},2),1:size(maps{1},1));
expX = expX.*binSize;       expY = expY.*binSize;          % Convert the x,y from bin counts to actual position
expX = expX-(binSize/2);    expY = expY-(binSize/2);       %  ..
expX = expX(visInd);        expY = expY(visInd);           % Convert to 1-D vector list of bins with occupancy>0, format (nBin,1)
% .. and we also need a linearised dwell time probability map, using the same 'visInd' %
posMap = posMap(visInd);
posMap = posMap ./ nansum(posMap(:));   % If we ever use the rate maps of one trial to decode the spikes of another, there may be NaNs in pos map (visited positions do not match up).



%------------------------------------ Some debugging/checks on position and spike data ---------------------------------------------------------------%
timeWinsInUse = sum(~isnan(actualX)) / nWin * 100;
timeWinsWithoutSpikes = sum(sum(actualFiring,2)==0 & ~isnan(actualX')) / nWin * 100;
% disp(['spk_positiondecode: ' num2str(size(maps,2)) ' cells, ' num2str(timeWinsInUse,'%3.2f') ' time wins in use, ' num2str(timeWinsWithoutSpikes,'%3.2f') ' of these without spikes.']);
% % % % Remove time windows with no spikes, if requested %
% % % if 0
% % %     actualFiring = actualFiring(~timeWinsWithoutSpikes,:);
% % %     actualX = actualX(~timeWinsWithoutSpikes);
% % %     actualY = actualY(~timeWinsWithoutSpikes);
% % %     nWin = sum(~timeWinsWithoutSpikes);
% % % end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Position decoding (sub-function call to specific method) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(decodeMethod,'MLE')  
    if ~usePriorPosition
        [pVals, pSpk] = MLEDecoding(expFiring,actualFiring,posMap,0); % 'pVals' is the key variable, the list of the probabilities of being in each pos bin, during each time window. Format is (1:nWin, 1:nBin).
    else
        pVals = MLEDecoding(expFiring,actualFiring,posMap,1,expX,expY,scaledAlpha,alphaVals,size(maps{1}),binSize,visInd);
    end   
elseif strcmp(decodeMethod,'directBasis')   
    pVals = directBasisDecoding(expX,expY,actualFiring);         % For direct basis, 'pVals' should actually be interpreted as the dot products of cell population firing rate vectors. Format is still (1:nWin, 1:nBin).
end
pVals = pVals ./ repmat(  nansum(pVals,2), 1, size(pVals,2)  );  % For all types of estimate, the pVals should be scaled such that, in each time window, the p for all spatial bins sums to 1.
pSpk = pSpk ./ repmat(  nansum(pSpk,2), 1, size(pSpk,2)  );  %

% For dir data, we can also smooth the p vals before caculating the best position %
% (Not done this yet for 2D data, the way I have collapsed xy data to a row vector makes this a headache %
if smoothPMode
    k = ones(1,varargin{2}) ./ varargin{2};
    pVals = imfilter(pVals,k,'circular');
end


% ----------------------------------------- Calculate distance error in decoding---------------------------------------------------------------------%
[maxPVals,mostLikelyBins] = nanmax(pVals,[],2);    % This tells you which is the most likely bin, and how likely you were to be there.
                                                   % Further note: if we ever use the rate maps of one trial to decode the spikes of another, there may be NaNs in pVals (as visited positions do not match up).
decodeErrors = sqrt(  (actualX' - expX(mostLikelyBins)).^2   +   (actualY' - expY(mostLikelyBins)).^2   );    % This still has the format (1:nWin)
if isDirMode
    % For dir data, errors so far will be in the range 0-360. Errors > 180 need to be mapped onto range 0-180. 
    decodeErrors( decodeErrors>180 )  = abs(  decodeErrors( decodeErrors>180 ) - 360  );   % -360 gets the circular equivalent, this is in range -180-0, but we only want the absolute error.
end


% --------- Recover p(pos|spk) at each actual visited position (i.e. never mind the max pVal, at the what was the pVal at the *actual* position of the rat, in each window T?) --------%
binnedY = ceil(actualY./binSize);                          %
binnedX = ceil(actualX./binSize);                          % All these lines here are to prevent rounding errors arising from CEIL at immediately following line.
binnedY(binnedY>size(maps{1},2)) = size(maps{1},2);        %
binnedX(binnedX>size(maps{1},1)) = size(maps{1},1);        %
actualPosInd = sub2ind(       size(maps{1}),    binnedY,     binnedX        );    % First, convert actualX and Y (format (1,nWin)) to a linear index into the rate map.
% Now, it gets complicated. We need an index that says where the rat was at each T, interms of the linear bin index 1:nVisBin,.
% Because of the way I have used visInd as the linear list of bins, rather than just a linearisation of the rate map,
% the indexing from the "cartesian" linear index in 'actualPosInd', to the one used in visInd requires the construction of an 'inverse vis ind', i.e.
% at each point in the "cartesian" linear index, where are we in 'visInd'?
invVisInd = zeros(numel(maps{1}),1);   % Pre-allocate
count=1;                                                            % Step through invVisInd (length=numel(map)), if bin is visisted 
for ii=1:length(invVisInd)                                          % (visInd=true), bump the count.
    if visInd(ii);  invVisInd(ii)=count;  count=count+1;  end       % After this, invVisInd is length=numel(map), but max=sum(visInd).
end                                                                 % When indexed by actualPosInd, it will map this, a straighforward linearisation of the rate map, onto visInd.
actualPosIsNan = isnan(actualPosInd);   % Because actualPosInd is an index, it can't have NaNs. (NaN=pos filtered due to edge or speed). Take them out here,
actualPosInd(actualPosIsNan) = 1;       % as format of actualVisPosInd doesn't change (always 1:nVisBin), we can put them back later.
actualVisPosInd = invVisInd(actualPosInd);       % Here, we actually do the conversion from "cartesian" linear index, to "visted" linear index.
actualVisPosIsZero = actualVisPosInd==0;   % A zero in 'actualVisPosInd' means that the mean actual postion in current T is a postion that is NOT actually visited (hence is not in rate 
actualVisPosInd(actualVisPosIsZero) = 1;   % map, and is zero in 'invVisInd'. We will write off these T, set to 1 for now, then set to NaN below.
actualVisPosIndLinear = sub2ind(  size(pVals), 1:size(pVals,1),  actualVisPosInd');
pValsActualPos = pVals(    actualVisPosIndLinear   );
pValsActualPos(actualPosIsNan | actualVisPosIsZero') = NaN;




%---------------------------------------------------Put output in nice structure---------------------------------------------------------------------%
% Positions %
rtn.errors = decodeErrors;
rtn.x = actualX';   % Flip these round so that everything is a column vector
rtn.y = actualY';   %
rtn.xPred = expX(mostLikelyBins);
rtn.yPred = expY(mostLikelyBins);
% Probabilities %
rtn.pValsBest = maxPVals;
rtn.pValsActualPos = pValsActualPos;
rtn.pVals = pVals;
rtn.pSpk = pSpk;
% Some parameters of the decode %
rtn.visInd = visInd;
rtn.speed = speedByWin';   % Flip these round so that everything is a column vector
rtn.perTimeWinInUse = propNanPerT';
rtn.timeWinsInUse = timeWinsInUse;
rtn.timeWinsWithoutSpikes = timeWinsWithoutSpikes;
rtn.timeWinsWithoutSpikesInd = nansum(actualFiring,2)==0;


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
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pVals pSpk] = MLEDecoding(expFiring,actualFiring,posMap,usePriorPosition,expX,expY,scaledAlpha,alphaVals,mapSize,binSize,visInd)
%----------------------------- Calculate most likely positions (basic MLE method, no prior) ---------------------------------------------------- ----%
% First, we want to vectorise, so create matching 3D matrices for expected and actual firing %
nWin = size(actualFiring,1);
temp = shiftdim(expFiring',-1);
expFiring3D = repmat(temp,[nWin 1 1]);
actualFiring3D = repmat(actualFiring,[1 1 size(expFiring,1)]);   % Format is (actual,nCell,Expected), i.e. (1:nWin, 1:nCell, 1:nBin), for both.
% Now calculate probabilities %
pValsEachCell    = (   (expFiring3D.^actualFiring3D)  ./  factorial(actualFiring3D)   ) .* exp( - expFiring3D );  
pVals = prod( pValsEachCell, 2 );   % This variable defines P(spk|pos), assuming independent poisson spiking for each cell. See Zhang et al 1998, Eq 35. Format is (1:nWin, 1:nBin).
pVals = squeeze( pVals );           % 'pVals' is the key variable, the list of the probabilities of being in each pos bin, during each time window.
pSpk = pVals;   % Keep a record of this, which is P(spk|pos) (i.e. what is the probability of the actual spiking seen in in T, given the rate map).
% We now multiply by the overall probability of being in any given position (i.e. the pos map)
% pVals then becomes P(spk|pos)*P(pos), which is equivalent to P(pos|spk), by Bayes rule. (If you ignore P(spk), which we will, by normalising P(pos|spk) such that it sums to 1 across the environment, see Zhang et al p)
pVals = pVals .* repmat(posMap',size(actualFiring,1),1);
% disp('NOT BAYESIAN - MLE');



%------------------- Combine 'online' position estimate with prior probability based on previous position -----------------------%
if usePriorPosition
    % Now need to cycle through time windows one-by-one, getting estimate on the basis of current ML and prior %
    pTempCurr = zeros(mapSize);
    pTempPrior = zeros(mapSize);
    pValsPrior = zeros(size(pVals));    % To store all prior P dists, for debugging purposes rather than for calculation.
    [X,Y] = meshgrid(1:mapSize(2),1:mapSize(1));   % Grid for construction of 2-D gaussian .. 
    X = X.*binSize;      Y = Y.*binSize;          %  .. convert the x,y from bin counts to actual position (in pix)
    X = X-(binSize/2);   Y = Y-(binSize/2);       %  ..
    for ii=1:nWin
        pTempCurr(:)=0;  pTempPrior(:)=0;             % Reset prob dist variables to all zero.
        pTempCurr(visInd) = pVals(ii,:);              % Get the ML (spike-based) p distribution for the current time window. At this point, format changes from 1:nBin vector to 2D spatial map. 
        pTempCurr = pTempCurr./nansum(pTempCurr(:));  % pVals need to be scaled to become an actual prob distribution, i.e sum to 1.
        if ii==1
            pTempPrior=ones(size(pTempPrior));     % In first time bin, no info on prior, so flat distribution.
        else
            prevX=expX(predPosTemp);   prevY=expY(predPosTemp);   % Get the previous estimated position.
            if scaledAlpha;   
                alpha=alphaVals(ii-1);   
            else
                alpha=alphaVals;
            end
            exponent = (    ((X-prevX).^2)  +  ((Y-prevY).^2)    )    ./   (2*alpha*alpha); % 2D gaussian centred on previous estimated postion, (optional: with
            pTempPrior = (exp(-exponent));                                                  % standard deviation scaled by speed in previous time window).               
        end
        pTempPrior = pTempPrior ./ nansum(pTempPrior(:));   % scaled to become an actual prob distribution that sums to 1 over the whole box.
        pCurrTimesPrior = pTempCurr.*pTempPrior;
        pCurrTimesPrior = pCurrTimesPrior ./ nansum(pCurrTimesPrior(:));
        pVals(ii,:) = pCurrTimesPrior(visInd);             % New predicted position, incorporating prior prediction, is stored in variable 'pVals
        pValsPrior(ii,:) = pTempPrior(visInd);   % Store all prior P dists, for debugging purposes rather than for calculation.
        [~,predPosTemp] = max(pVals(ii,:));      % This is the current predicted position, will be needed to set the prior in next time window.
    end
end
                                            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corrVals = directBasisDecoding(expFiring,actualFiring)
% First, we want to vectorise, so create matching 3D matrices for expected and actual firing %
temp = shiftdim(expFiring',-1);
expFiring3D = repmat(temp,[nWin 1 1]);
actualFiring3D = repmat(actualFiring,[1 1 size(expFiring,1)]);   % Format is (actual,nCell,Expected), i.e. (1:nWin, 1:nCell, 1:nBin), for both.
corrVals = dot(expFiring3D,actualFiring3D,2);        % Take the dot product along the 1:nCell dimension (i.e. the dot product of the population vectors at each time window and spatial bin.
corrVals = squeeze(corrVals);     % Final format is (1:nWin, 1:nBin)








