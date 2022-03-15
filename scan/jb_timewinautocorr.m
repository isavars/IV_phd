function []=jb_timewinautocorr(varargin)
% Time-windowed spike triggered spatial (directional) auto-correlograms

prms.mode = 'auto';
prms.bin = 36;              % Binsize, in camera pixels for place, in degrees for direction.
prms.boxcar = 2;            % Use fspecial('gaussian', etc ) here for gaussian smooth, single integer A for AxA flat boxcar.
prms.dwellThresh = 5;          % AC overlap threshold, units = sec
prms.T = 5:2:35;   % Time windows
prms.preMadeData = [];
prms.turnLimit = [];
for ii=1:2:length(varargin)-1
    prms.(varargin{ii}) = varargin{ii+1};
end

SD=gss;
data=evalin('base',SD.selData{1});
trInd=SD.selTrial;
cellInd=SD.selCell;

dirMapsWholeTrial = rates_main(data, rates_params('space','dir'));

if strcmp(prms.mode,'auto')
    % Autocorrs %
    hFig = gra_multiplot(length(cellInd),length(trInd)*2);   hAxArray=getappdata(hFig,'axesHandles');   figure(hFig);   colormap(gra_tintcolormap(50));
    for ii=1:length(trInd)
        trTemp=data.trials(trInd(ii));
        for jj=1:length(cellInd)
            % Plot regular Dir map for comparison %
            gra_plotmap(dirMapsWholeTrial{trInd(ii),cellInd(jj)},'parent',hAxArray( jj, ii*2) );
            % Do not create spike triggered maps if no spikes %
            if isempty(trTemp.cells(cellInd(jj)).st);    continue;      end
            % Run AC and PLot %
            if isempty(prms.preMadeData)
                % Run time-win corr on raw data %
                [tWinDirMaps,binVect] = spk_spktrigtimewincorr_dir(trTemp.dir, trTemp.cells(cellInd(jj)).st,trTemp.cells(cellInd(jj)).st, prms.T, trTemp.sample_rate, prms.bin, prms.boxcar, prms.dwellThresh);
                combineTimeLagsAndPlot(prms,hAxArray( jj,  (ii*2)-1),tWinDirMaps,binVect);
            else
                % Use pre-made corr from mega-array, made for all data %
                combineTimeLagsAndPlot(prms,hAxArray( jj,  (ii*2)-1),prms.preMadeData,SD.selData{1},trInd(ii),cellInd(jj), cellInd(jj));
            end
        end
    end
    gra_multilabel(hFig,'title',{SD.selData{1}, ['bin=' num2str(prms.bin)],  ['dwellThresh=' num2str(prms.dwellThresh)]});
    colLabels = cell(1,length(trInd)*2);   
    colLabels(1,2:2:end) = getdatanames(data,'trial','long',SD.selTrial);
    gra_multilabel(hFig,'col',colLabels);
    gra_multilabel(hFig,'row',getdatanames(data,'cell','short',SD.selCell));
    
elseif strcmp(prms.mode,'cross')
    % Cross-corrs %
    trTemp=data.trials(trInd(1));
    hFig = gra_multiplot(length(cellInd)+1,length(cellInd)+1);   hAxArray=getappdata(hFig,'axesHandles');   figure(hFig);  colormap(gra_tintcolormap(50));
    for ii=1:length(cellInd)
        for jj=1:length(cellInd)
            if jj==ii;   continue;    end
            if isempty(prms.preMadeData)
                % Run time-win corr on raw data %
                [tWinDirMaps,binVect] = spk_spktrigtimewincorr_dir(trTemp.dir, trTemp.cells(cellInd(ii)).st,trTemp.cells(cellInd(jj)).st, prms.T, trTemp.sample_rate, prms.bin, prms.boxcar, prms.dwellThresh);
                combineTimeLagsAndPlot(prms,hAxArray(jj+1,ii+1),tWinDirMaps,binVect);
            else
                % Use pre-made corr from mega-array, made for all data %
                combineTimeLagsAndPlot(prms,hAxArray(jj+1,ii+1),SD.selData{1},trInd(1),cellInd(ii), cellInd(jj));
            end
        end
    end
    for ii=1:length(cellInd)
        gra_plotmap(dirMapsWholeTrial{trInd(1),cellInd(ii)},'parent',hAxArray( 1, ii+1) );
        gra_plotmap(dirMapsWholeTrial{trInd(1),cellInd(ii)},'parent',hAxArray( ii+1, 1) );
    end
    gra_multilabel(hFig,'row',[{''} getdatanames(data,'cell','short',SD.selCell)]);
    gra_multilabel(hFig,'col',[{''} getdatanames(data,'cell','short',SD.selCell)]);
    gra_multilabel(hFig,'title',{SD.selData{1},trTemp.trialname, ['bin=' num2str(prms.bin)],  ['dwellThresh=' num2str(prms.dwellThresh)]});
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=combineTimeLagsAndPlot(prms,hAx,varargin)   % tWinDirMaps,binVect)
% Combine separate time-lag maps into an image matrix %
% First parse input: has this function been supplied with just one corr, or rather the mega-array of pre-made corrs for the whole dataset? %
if length(varargin)==2
    % Just one corr supplied ..  %
    tWinDirMaps = varargin{1};   binVect=varargin{2};
    % Also need to process the raw output of SPK_SPKTRIGTIMEWINCORR (cell array, one cell for each time lag) to a 2D array, fomat (timeLag,spatBin) %
    tWinMapAll = nan(length(tWinDirMaps), length(binVect));                             
    for ii=1:length(tWinDirMaps)
        tWinMapAll(ii,:) = tWinDirMaps{ii}';
    end
else
    % Mega-array supplied - need to use suppled indices to extract the correct one %
    dsInd = strcmp(prms.preMadeData.datasetName,varargin{1});
    tWinMapAll = prms.preMadeData.corrs{ dsInd }{ varargin{2} }{ varargin{3}, varargin{4} };
    binVect = prms.preMadeData.binVects{ dsInd }{ varargin{2} }{ varargin{3}, varargin{4} };
    % As mega-array is calculated with no dwell threshold, need to apply this, using matching dwell map %
    posForCorr = prms.preMadeData.pos{ dsInd }{ varargin{2} }{ varargin{3}, varargin{4} };
    tWinMapAll(  posForCorr < prms.dwellThresh  ) = NaN;
end

if isempty(binVect);   return;   end

% Remove data from beyond a certain number of turns, if requested %
if ~isempty(prms.turnLimit)
    turnLimInd = repmat(   (abs(binVect)./360)>prms.turnLimit, size(tWinMapAll,1), 1)    &    ~isnan(tWinMapAll);
    tWinMapAll(  turnLimInd  ) = 0;
end

% Plot Time-Win map %
gra_plotmap(tWinMapAll,'parent',hAx,'text_pos','none');
axis(hAx,'normal','on');
if strcmp(prms.mode,'cross')
    tWinList = [fliplr(-prms.T) prms.T];
else
    tWinList = prms.T;
end
yTickInd = 1:5:length(tWinList);   set(hAx,'ytick',yTickInd);   set(hAx,'yticklabel',tWinList(yTickInd));
mapAllCenInd = find(  binVect==0  );    
nTurns = floor(  (mapAllCenInd*prms.bin)  / 360  );
if nTurns<2;   nTurnsAxLim=2;    else    nTurnsAxLim=nTurns+1;    end
set(hAx,'xlim',([-1 1].*(nTurnsAxLim).*(360/prms.bin))+mapAllCenInd);
hold(hAx,'on');
xTickInd = -nTurnsAxLim:nTurnsAxLim;  set(hAx,'xtick',mapAllCenInd - 0.5 + (xTickInd*360)./prms.bin);   set(hAx,'xticklabel',xTickInd);
for kk=2:length(xTickInd)-1
    plot(hAx, mapAllCenInd - 0.5 + [1 1].*((xTickInd(kk)*360)/prms.bin), get(hAx,'ylim'), 'k:');
end
     




















        
        
        