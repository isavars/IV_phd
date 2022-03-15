
function []=gra_spktrigtimewincorr(varargin)
% Plot spike-triggered time-windowed spatial correlograms
%
%       gra_spktrigtimewincorr()


% Input option defaults: for this function %
opt.mode = 'cross';
opt.T = 10;
opt.binsize = 10;
opt.nColour = 20;
opt.trimEdges = 1;
% Input option defaults: for SPK_SPKTRIGTIMEWINCORR %
opt.dwellThresh = 2;
opt.boxcar = 5;
opt.posShuffleActive = 1;
opt.posShuffleMinOffset = 20;
opt.nPosShuffles = 100;
opt.normaliseToShuf = 1;
opt.temporalDownSamp = 5;
opt.showAs1D = 0;               % Collapse x and y to 1D  radial displacement.
opt.offsetTimeWindow = 0;
opt.isDirData = 1;
% Parse caller-supplied input options %
for ii=1:2:length(varargin)-1
    opt.(varargin{ii}) = varargin{ii+1};
end

SD=gss;
data=evalin('base',SD.selData{1});
trInd=SD.selTrial;
cellInd=SD.selCell;
cMap = [1 1 1; jet(opt.nColour-1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main cell vs cell loop %
if strcmp( opt.mode, 'cross' )
    % Each cell against another witinh a trial %
    if length(trInd)>1; errordlg('When using ''cross'' mode, you can only select one trial at a time', 'Error in GRA_SPKTRIGTIMEWINCORR'); end
    
    hFig = gra_multiplot(length(cellInd),length(cellInd));   hAx=getappdata(hFig,'axesHandles');        colormap(hFig,cMap);
    trTemp = data.trials(trInd);
    hWait=waitbar(0,'Please wait');
    for ii=1:length(cellInd)
        for jj=1:length(cellInd)
            if opt.isDirData
                [tWinMap,binVect] = spk_spktrigtimewincorr(trTemp.dir,     [], trTemp.cells(cellInd(ii)).st, trTemp.cells(cellInd(jj)).st, opt.T, trTemp.sample_rate, opt.binsize, opt);
            else
                [tWinMap,binVect] = spk_spktrigtimewincorr(trTemp.x, trTemp.y, trTemp.cells(cellInd(ii)).st, trTemp.cells(cellInd(jj)).st, opt.T, trTemp.sample_rate, opt.binsize, opt);
            end
            plotCorr( tWinMap, binVect, hAx(ii,jj), opt );
            if ii==jj;
                set(hAx(ii,jj),'Visible','on','Box','on','xtick',[],'ytick',[]);   % Highlight the diagonal line of autocorrs by boxing the axes
            end
            waitbar(  (jj + ((ii-1)*length(cellInd))) / (length(cellInd)^2), hWait );
        end 
    end
    gra_multilabel(hFig,'row',getdatanames(data,'cell','short',cellInd));
    gra_multilabel(hFig,'col',getdatanames(data,'cell','short',cellInd));
    titleStr = {SD.selData{1},data.trials(trInd).trialname,num2str(opt.T)};

elseif strcmp( opt.mode, 'auto' )
    % Each cell against itself, for mulitple trials %
    hFig = gra_multiplot(length(cellInd),length(trInd));   hAx=getappdata(hFig,'axesHandles');         colormap(hFig,cMap);
    hWait=waitbar(0,'Please wait');
    for ii=1:length(trInd)
        for jj=1:length(cellInd)
            trTemp = data.trials(trInd(ii));
            if opt.isDirData
                [tWinMap,binVect] = spk_spktrigtimewincorr(trTemp.dir,  []  ,trTemp.cells(cellInd(jj)).st,trTemp.cells(cellInd(jj)).st, opt.T, trTemp.sample_rate, opt.binsize, opt);
            else
                [tWinMap,binVect] = spk_spktrigtimewincorr(trTemp.x,trTemp.y,trTemp.cells(cellInd(jj)).st,trTemp.cells(cellInd(jj)).st, opt.T, trTemp.sample_rate, opt.binsize, opt);
            end
            plotCorr( tWinMap, binVect, hAx(jj,ii), opt );
            waitbar(  (jj + ((ii-1)*length(trInd))) / (length(cellInd)*length(trInd)), hWait );
        end
    end
    gra_multilabel(hFig,'row',getdatanames(data,'cell','short',cellInd));
    gra_multilabel(hFig,'col',getdatanames(data,'trial','long',trInd));
    titleStr = {SD.selData{1},num2str(opt.T)};
end
gra_multilabel(hFig,'title',titleStr);
close(hWait);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=plotCorr(map,binVect,hAx,opt)
% Sub-function to plot corr maps %
if min(size(map))==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1D maps: Directional data and 2D maps collapsd to radial distance %
    plot(hAx,binVect,map,'b-');
    hold(hAx,'on');
    plot(hAx,get(hAx,'xLim'),[0 0],'k-');   
    if opt.normaliseToShuf && opt.posShuffleActive
        plot(hAx,get(hAx,'xLim'),[-2 -2],'k:');   plot(hAx,get(hAx,'xLim'),[2 2],'k:');
    end
    % Set axes limits %
    mapLim = max( abs([nanmin(map) nanmax(map)]) );
    if mapLim<2.5 && opt.normaliseToShuf && opt.posShuffleActive;  mapLim=2.5;  end
    set(hAx,'yLim',[-mapLim, mapLim]);
    if opt.isDirData
       set(hAx,'xlim',[-360 360],'xtick',-360:90:360,'xgrid','on'); 
    end
    
%     hold(hAx(ii,jj),'on');
%     dirsWithData = binVect{jj}(   tWinMap{jj}>0   );
%     if isempty(dirsWithData) || max(dirsWithData)<450    % Is empty when there are no spikes, no spikes in window, or no pos>dwellThresh in map.
%         set(hAx(ii,jj),'xlim', [-450 450]);
%     else
%         set(hAx(ii,jj),'xlim', [min( dirsWithData ) max(dirsWithData)]);
%     end
%     plot(hAx(ii,jj),[-360 -360],get(hAx(ii,jj),'ylim'),'k:');
%     plot(hAx(ii,jj),[360 360],get(hAx(ii,jj),'ylim'),'k:');
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2D maps                                                           %
    % If requested, remove bits of the image that are not connected to the centre %
    if opt.trimEdges
        labels = bwlabel( ~isnan(map) , 4 );
        cenLabel = labels(round(numel(labels)/2));
        map( labels~=cenLabel ) = nan;
    end
    % Draw image %
    minRate=nanmin(map(:));  % Get the min and max of the real data 
    maxRate=nanmax(map(:));  %
    map(isnan(map)) = minRate - ((maxRate-minRate)/(opt.nColour-1.5));  % Set the nans (unvisited) to real values, designed to be the lowest colour step. This has been set as [1 1 1].
    imagesc(map,'parent',hAx);
    axis(hAx,'equal','off');
    % Print min and max rate %
    text('string',num2str(maxRate,'%1.1f'),'units','normalized','position',[0 1],'HorizontalAlignment','left','FontUnits','normalized','VerticalAlignment','cap','fontsize',0.15,'parent',hAx);
    text('string',num2str(minRate,'%1.1f'),'units','normalized','position',[1 1],'HorizontalAlignment','right','FontUnits','normalized','VerticalAlignment','cap','fontsize',0.15,'parent',hAx);
    % Cross-hairs through the centre point %
    hold(hAx,'on');
    plot(hAx,[1 1].*find(binVect==0), get(hAx,'yLim'), 'k-');
    plot(hAx,get(hAx,'xLim'), [1 1].*find(binVect==0), 'k-');
end







