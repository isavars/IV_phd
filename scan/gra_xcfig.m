function []=gra_xcfig(SD,varargin)
% Multiple plot of temporal cross- or auto-correlograms:
%   gra_xcfig(SD)
%
% To adjust the Correlogram window or bin, add the arguments:
%       gra_xcfig(SD,'xcBin', bin, 'xcWin', win);
%
% To plot an array of cross-correlograms, select one trial only, and call:
%   gra_ac_fig(SD,'corrMode','cross')
%
% To add the spatial maps to a cross-corr plot, add the arguments ..
%    gra_acfig( ... 'addSpatialMap',1)
%

% opengl SOFTWARE     % To prevent hardware bug on Dell T1600s, which flips text upside down.

data = evalin('base', SD.selData{1});
trInd=SD.selTrial;
cellInd=SD.selCell;
if length(trInd)>1;  disp('GRA_XCFIG: need to select one trial at a time.');    trInd=trInd(1);  end

% Parse input parameters %
prms.corrMode='auto';
prms.xcBin=0.001;
prms.xcWin=0.1;
prms.normaliseToShuf = 0;
prms.normaliseToShufMode = 'norm';
prms.nShuffles = 100;
prms.shuffleMinOffset = 20;
prms.addSpatialMap=0;
prms.plotSigThr = [1];   %  If not empty, plots horizontal lines at the specified points. Most useful when used with spike-shift normalised XCs.
prms.showMeanCenPart = 0.5;
prms.spikeProb = 0;  % With this option, Y-Axis will show the probability of a reference cell spike being followed by a test cell spike. (xcorr normalisation is 'none', and xc is divided by N reference spikes), e.g. Dupret et al 2013. 
for ii=1:2:length(varargin)
    prms.(varargin{ii}) = varargin{ii+1};
end

%%%%%%%% REMOVED JULY 2017. TW - all these old-school ways of normalising xc no longer necessary after writing proper multi-time shift code, which is in SPK_CROSSCORR
% prms.xCorrSubtraction='none'  ; % 'multiShift'; % 'shiftedSpikeTrain'; %  % '30-50ms'; %   % 'shiftedSpikeTrain';      % Method for subtractive normalisation of xcorrs - see Csicivari et al  1998 or Dupret et al 2013.


% Run the correlograms %
switch prms.corrMode
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'auto'
        hFig = gra_multiplot(length(cellInd),length(trInd),'graphborder',repmat(0.2, [1 4]));   
        ud=get(hFig, 'userdata');   ax=ud.axes_handle_array;
        for ii=1:length(trInd)
            for jj=1:length(cellInd)
                st=data.trials(trInd(ii)).cells(cellInd(jj)).st;
                prms.plot = ax(jj,ii);
                spk_crosscorr(st, 'AC', prms.xcBin, prms.xcWin, data.trials(trInd(ii)).dur, prms);
                set(ax(jj,ii),'yticklabel','');
            end
        end
        gra_multilabel(hFig, 'row', getdatanames(data,'cell','short',cellInd));
        gra_multilabel(hFig, 'col', getdatanames(data,'trial','short',trInd));
        
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'cross'
        if ~prms.addSpatialMap
            hFig = gra_multiplot(length(cellInd),length(cellInd),'axisborder',repmat(0.3, [1 4]),'figborder',[0.8 0.5 0.5 0.5]);
        else
            hFig = gra_multiplot(length(cellInd)+1,length(cellInd)+1,'axisborder',repmat(0.3, [1 4]),'figborder',[0.8 0.5 0.5 0.5]);            
        end
        ud=get(hFig, 'userdata');   axArray=ud.axes_handle_array;
        for ii=1:length(cellInd)
            for jj=1:length(cellInd)
                if ii==jj
                    axis(axArray(ii,jj),'off');
                    continue
                end
                ax=axArray(ii,jj);
                stA=data.trials(trInd).cells(cellInd(jj)).st;
                stB=data.trials(trInd).cells(cellInd(ii)).st;
                prms.plot = ax;
                [xc,lags] = spk_crosscorr(stB,stA,prms.xcBin,prms.xcWin,data.trials(trInd).dur,prms);      % Note that if the 'reference' spike train goes as second argument, then the intuitively meaningful lag values come out as the second half of the xcorr.

%%%%%%%% REMOVED JULY 2017. TW - all these old-school ways of normalising xc no longer necessary after writing proper multi-time shift code, which is in SPK_CROSSCORR
%                 % ---------------- Subtractive normalisation --------------------- %
%                 if strcmp(prms.xCorrSubtraction,'shiftedSpikeTrain')
%                     stB = stB+0.1;
%                     acShift=runXCorr(stB,stA,prms.xcBin,prms.acWin,data.trials(trInd).dur,prms);
%                     ac = ac - acShift;
%                 elseif strcmp(prms.xCorrSubtraction,'multiShift')
%                     nShift=10;
%                     acShift = nan(nShift,length(ac));
%                     for kk=1:nShift
%                         stB = stB+kk;
%                         acShift(kk,:)=runXCorr(stB,stA,prms.xcBin,prms.acWin,data.trials(trInd).dur);
%                     end
%                     ac = ac - nanmean(acShift)';
%                 elseif strcmp(prms.xCorrSubtraction,'30-50ms')
%                     ac = ac - mean(ac(30:50));  % ATTENTION! Hard coded time values here assume xc bin=1ms, xc win >= 50ms. !!!!
%                 end
%                 % ------ Are we making a 'probablility of following spike' plot? ------- %
%                 if prms.spikingProbYAxis
%                     if prms.xcBin~=1;   disp('Warning: you cannot trust the Y-axis (Spike Prob) of these plots (as xcorr time bin >1ms)');   end
%                     ac = ac ./ length(stA);
%                 end
%                 % ----------------------- Plot ------------------------------------ %
%                 plotXCorr(ac,prms.xcBin,prms.acWin,ax);
%%%%%%%% REMOVED JULY 2017

                % --------------- Significance threshold -------------------------- %
                if ~isempty( prms.plotSigThr )
                    hold(ax,'on');
                    for kk=1:length(prms.plotSigThr)
                        plot(ax,get(ax,'xlim'),[1 1].*prms.plotSigThr(kk),'m:');
                    end
                end
                % Mean of central part - draw lines to show, and calculate and print mean on plot %
                if ~isempty( prms.showMeanCenPart )
                    t = prms.showMeanCenPart;      rep=[-1 1];
                    for kk=1:2;   plot(ax,[t t].*rep(kk),get(ax,'ylim'),'b:');     end
                    meanCenPart = nanmean( xc( (lags*prms.xcBin)>=-prms.showMeanCenPart & (lags*prms.xcBin)<=prms.showMeanCenPart ) );
                    text('units', 'normalized', 'position', [0 1], 'HorizontalAlignment', 'left', 'string', num2str(meanCenPart,'%3.2f'), 'FontUnits', 'normalized', 'VerticalAlignment', 'cap', 'fontsize', 0.1,'parent',ax);
                end
                if ~( prms.spikeProb ||  prms.normaliseToShuf )
                    set(ax,'yticklabel','');
                end
            end
        end
        
        if prms.addSpatialMap
            for ii=1:length(cellInd)
                gra_plotmap(data.rate_maps( SD.selMap ).maps{trInd,cellInd(ii)},'handle',axArray(end,ii));
                gra_plotmap(data.rate_maps( SD.selMap ).maps{trInd,cellInd(ii)},'handle',axArray(ii,end));
            end
        end
        
        str = getdatanames(data,'cell','short',cellInd);
        for ii=1:length(str);  str{ii} = [str(ii) {'(Ref)'}];  end 
        gra_multilabel(hFig, 'col', str);
        str = getdatanames(data,'cell','short',cellInd);
        for ii=1:length(str);  str{ii} = [str(ii) {'(Test)'}];  end 
        gra_multilabel(hFig, 'row', str);
        gra_multilabel(hFig, 'title', {data.trials(trInd).trialname, num2str(prms.xcWin), num2str(prms.xcBin)});
        
end

% opengl AUTOSELECT

%%%%%%%% REMOVED JULY 2017. TW - why did this function have its own plotting and cross-corr code? Now calls SPK_CROSSCORR.
% %-------------------------------------------------------------------------%
% function [ac]=runXCorr(stA,stB,xcBin,acWin,trDur,prms)
% %%% Do auto-correlation and plot %%%
% % Make spike train vector %
% % spkTrHistInd=ceil((st.*1000)/xcBin);     % As spkTr in sec, xcBin in ms
% % spkTrHist=repmat(0,1,ceil((trDur.*1000)./xcBin));
% % spkTrHist(spkTrHistInd)=1;
% spkTrHistA=histc(stA,0:(acBin/1000):trDur);
% spkTrHistB=histc(stB,0:(acBin/1000):trDur);
% % AC %
% [ac, lags]=xcorr(spkTrHistA,spkTrHistB,(acWin/acBin),'none');
% % Select only the relevant short window %
% if strcmp( prms.corrMode, 'auto' );   ac(lags==0) = NaN;   end
% winInd= lags>=-(acWin/acBin) & lags<=(acWin/acBin);
% ac=ac(winInd);
% 
% %------------------------------------------------------------------------%
% function [] = plotXCorr(ac,acBin,acWin,hAxis)
% %plot(hAxis,acBin:acBin:acWin,ac,'k-');
% bar(hAxis,-acWin:acBin:acWin,ac,'histc');
% shading(hAxis,'flat');
% if max(ac)==0; ylim = [0 1];  else   ylim = [min(ac) max(ac)];   end
% set(hAxis,'xlim',[-acWin acWin],'ylim',ylim,'xtick',[-acWin round(-acWin/2) 0 round(acWin/2) acWin]);
% set(hAxis,'tickdir','out');



