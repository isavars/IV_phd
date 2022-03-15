function [rtn] = spk_MLEpositiondecode(T,maps,spikeTimes,xSamp,ySamp,trialDur,posSampRate,binSize,usePriorPosition,varargin)

% HELP
%       [rtn] = spk_MLEpositiondecode(T,maps,spikeTimes,xSamp,ySamp,trialDur,posSampRate,binSize,usePriorPositionFlag)
%       [rtn] = spk_MLEpositiondecode(T,maps,spikeTimes,xSamp,ySamp,trialDur,posSampRate,binSize,usePriorPositionFlag,alpha)
%       [rtn] = spk_MLEpositiondecode(T,maps,spikeTimes,xSamp,ySamp,trialDur,posSampRate,binSize,usePriorPositionFlag,speed,alphaSpeedScaling)
%
% Input 'binSize' is the number of pixels in each bin, assuming that you have passed the x and the y in pixels as units.
%
% rtn.errors = decodeErrors;
% rtn.probs = maxPVals;
% rtn.x = actualX;
% rtn.y = actualY;



%------------------------------------------------ Parse Input ------------------------------------------------------------------%
% Get the alpha for the prior 2D gaussian. Can either be fixed, or scale with speed, depening on user input %
scaledAlpha=0;
if usePriorPosition
    if length(varargin)==1
        scaledAlpha = 0;
        alpha = varargin{1};    % Fixed
    elseif length(varargin)==2
        scaledAlpha = 1;
        speed = varargin{1};
        alphaSpeedScale = varargin{2};
    end
end
debug=0;
% disp(['Use Prior Position = ' num2str(usePriorPosition)]);



%-------------------------------- Mean the position and speed data into time windows -------------------------------------------%
% If sample rate=50, simple reshape. If DACQ1 (46.875Hz), a bit more long-winded %
nWin = floor(trialDur/T);
if posSampRate==50 && rem(T,1/posSampRate)==0
    % In the case when the sample rate is 50, and the time window is divisible by 1/50, simple reshape required, but may need to cut off 'overhanging' data that doesn't make a full time window % 
    xSamp=xSamp(1:nWin*T*posSampRate);
    ySamp=ySamp(1:nWin*T*posSampRate);
    actualX = nanmean(   reshape(xSamp,posSampRate*T,nWin)    ,  1  );
    actualY = nanmean(   reshape(ySamp,posSampRate*T,nWin)    ,  1  );
    if scaledAlpha;
        speed=speed(1:nWin*T*posSampRate);
        speedByWin = nanmean(   reshape(speed,posSampRate*T,nWin)    ,  1  );   
    end    
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
        if scaledAlpha;   speedTemp(1:sum(ind),ii) = speed(ind);   end
    end
    actualX = nanmean(xTemp);
    actualY = nanmean(yTemp);
    if scaledAlpha;   speedByWin = nanmean(speedTemp);   end    
end
if scaledAlpha;   speedForAlpha = (speedByWin .* alphaSpeedScale(2)) + alphaSpeedScale(1);   end  % Actual speeds need to be linearly re-scaled into a range of alpha values. See Zhang et al 1998.



%--------------------------------------------------------------------------------------------------------------------------------%
% Generate 2D matrix of actual observed firing during each time window, format (1:nWin,1:nCell), where nWin = trial_dur/T %
actualFiring = nan(nWin,length(maps));  % Preallocate
for ii=1:length(spikeTimes)
    temp=histc(spikeTimes{ii},0:T:trialDur);
    actualFiring(:,ii) = temp(1:end-1);
end



%--------------------------------------------------------------------------------------------------------------------------------%
%%% Generate 2D matrix of expected firing values, with format (1:nBin,1:nCell) %%%
visInd = ~isnan(maps{1});
expFiring = nan(     sum(visInd(:)),    length(maps)   );     % Preallocate - N Rows is number of visited bins, rather possible total in env.
% Loop through supplied maps, linearising each map into a column of the expected firing array %
for ii=1:length(maps)
    expFiring(:,ii) = maps{ii}(visInd);        % Lineraise using visInd.
    expFiring(:,ii) = expFiring(:,ii) .* T;    % Multiply by T to convert from Hz to expected spikes in T (time window)
end
% We also need an x and y co-ord lookup, to be able to retrive error in 2D %
[expX expY] = meshgrid(1:size(maps{1},2),1:size(maps{1},1));
expX = expX.*binSize;       expY = expY.*binSize;          % Convert the x,y from bin counts to actual position
expX = expX-(binSize/2);    expY = expY-(binSize/2);       %  ..
expX = expX(visInd);        expY = expY(visInd);





%--------------------------------------------------------------------------------------------------------------------------------%
% Some debugging/checks on position and spike data %
% if sum(isnan(actualX))>0 || sum(isnan(actualY))>0
%     dbstop in spk_MLEpositiondecode at 104
% end
% Check how many time windows have no spikes at all, and report this (as these time bins cause problems,  %
% the expected position always 'defaults' to the position bin with the lowest overall expected firing) %
timeWinsWithoutSpikes = sum(actualFiring,2)==0;
disp(['spk_MLEpositiondecode: ' num2str(size(maps,2)) ' cells, ' num2str(sum(timeWinsWithoutSpikes)) '/' num2str(nWin) '(' num2str(round((sum(timeWinsWithoutSpikes)/nWin)*100))  '%) time windows with no spikes.']);
% % % % Remove time windows with no spikes, if requested %
% % % if 0
% % %     actualFiring = actualFiring(~timeWinsWithoutSpikes,:);
% % %     actualX = actualX(~timeWinsWithoutSpikes);
% % %     actualY = actualY(~timeWinsWithoutSpikes);
% % %     nWin = sum(~timeWinsWithoutSpikes);
% % % end




%--------------------------------------------------------------------------------------------------------------------------------%
%%% Calculate most likely positions %%%
% First, we want to vectorise, so create matching 3D matrices for expected and actual firing %
temp = shiftdim(expFiring',-1);
expFiring3D = repmat(temp,[nWin 1 1]);
actualFiring3D = repmat(actualFiring,[1 1 size(expFiring,1)]);   % Format is (actual,nCell,Expected), i.e. (1:nWin, 1:nCell, 1:nBin), for both.
% Now calculate probabilities %
%    pval_contrib    = ((tuning.^currK)./fact_currK) .* expon;
pValsEachCell    = (   (expFiring3D.^actualFiring3D)  ./  factorial(actualFiring3D)   ) .* exp( - expFiring3D );
pVals = prod( pValsEachCell, 2 );
pVals = squeeze( pVals );           % This is the key variable, the list of the probabilities of being in each pos bin, during each time window. Format is (1:nWin, 1:nBin).



%--------------------------------------------------------------------------------------------------------------------------------%
%%% Caswell suggested an extra step: smooth the probability distribtuions before taking         %%%
%%% most likely position. Computationaly pain in the aras, would need to reconfigure the nBin   %%%
%%% list for each time window into a 2D map, and smooth it (or is there a clever matlab way?)   %%%     


%------------------- Combine 'online' position estimate with prior probability based on previous position -----------------------%
if usePriorPosition
    % Now need to cycle through time windows one-by-one, getting estimate on the basis of current ML and prior %
    pTempCurr = zeros(size(maps{1}));
    pTempPrior = zeros(size(maps{1}));
    pValsPrior = zeros(size(pVals));    % To store all prior P dists, for debugging purposes rather than for calculation.
    [X,Y] = meshgrid(1:size(maps{1},2), 1:size(maps{1},2));   % Grid for construction of 2-D gaussian .. 
    X = X.*binSize;      Y = Y.*binSize;          %  .. convert the x,y from bin counts to actual position (in pix)
    X = X-(binSize/2);   Y = Y-(binSize/2);       %  ..
    for ii=1:nWin
        pTempCurr(:)=0;  pTempPrior(:)=0;             % Reset prob dist variables to all zero.
        pTempCurr(visInd) = pVals(ii,:);              % Get the ML (spike-based) p distribution for the current time window. At this point, format changes from 1:nBin vector to 2D spatial map. 
        pTempCurr = pTempCurr./nansum(pTempCurr(:));  % pVals need to be scaled to become an actual prob distribution, i.e sum to 1. Shouldn't have NaNs in, but I'm finding a bug with occasional NaNs in visited positions - possibly in adaptive smoothing?
        if ii==1
            pTempPrior=ones(size(pTempPrior));     % In first time bin, no info on prior, so flat distribution.
        else
            prevX=expX(predPosTemp);   prevY=expY(predPosTemp);   % Get the previous estimated position.
            if scaledAlpha;   alpha=speedForAlpha(ii-1);    end
            exponent = (    ((X-prevX).^2)  +  ((Y-prevY).^2)    )    ./   (2*alpha*alpha); % 2D gaussian centred on previous estimated postion, (optional: with
            pTempPrior = (exp(-exponent));                                                  % standard deviation scaled by speed in previous time window).               
        end
        pTempPrior = pTempPrior ./ nansum(pTempPrior(:));   % scaled to become an actual prob distribution, i.e sum to 1. Shouldn't have NaNs in, but I'm finding a bug with occasional NaNs in visited positions - possibly in adaptive smoothing?
        pTemp = pTempCurr.*pTempPrior;
        pTemp = pTemp ./ nansum(pTemp(:));
        pVals(ii,:) = pTemp(visInd);             % New predicted position, incorporating prior prediction, is stored in variable 'pVals
        pValsPrior(ii,:) = pTempPrior(visInd);   % Store all prior P dists, for debugging purposes rather than for calculation.
        [~,predPosTemp] = max(pVals(ii,:));      % This is the current predicted position, will be needed to set the prior in next time window.
    end
end
              
                               

%--------------------------------------------------------------------------------------------------------------------------------%
%%% Calculate distance error in decoding %%%
[maxPVals,mostLikelyBins] = max(pVals,[],2);    % This tells you which is the most likely bin, and how likely you were to be there.
decodeErrors = sqrt(  (actualX' - expX(mostLikelyBins)).^2   +   (actualY' - expY(mostLikelyBins)).^2   );    % This still has the format (1:nWin)




%----------------------------- Debug (plot probs, real and extimated positions) -------------------------------------------------%
if debug
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
    


%--------------------------------------------------------------------------------------------------------------------------------%
%%% Put output in nice structure %%%
rtn.errors = decodeErrors;
rtn.predProbs = maxPVals;
rtn.x = actualX;
rtn.y = actualY;
rtn.xPred = expX(mostLikelyBins);
rtn.yPred = expY(mostLikelyBins);
rtn.pVals = pVals;
rtn.visInd = visInd;
rtn.timeWinsWithoutSpikes = sum(timeWinsWithoutSpikes) / length(timeWinsWithoutSpikes);









