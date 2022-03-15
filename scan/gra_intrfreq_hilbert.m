function []=gra_intrfreq_hilbert(SD)
% Plots results of two related analyses:
%  1) Phases of spikes (derived from instantaneous phase, ie. hilbert transform).
%  2) Intrinsic frequency, based inter-theta cycle burst delays.

%%%%%%%%%%%%%%%%%%%%%%%%
bandWidth=3;            % This is actually half-band width (pass band is peak+-bandWidth)
speedLimit=2.5;
powerThr=[5 100];      % Inclusion cut-offs for power per cycle (%tiles of powerPerCyle for whole trial). Applied after speed filtering.
maxElapsedCycles=1;     % How many cycles to count across (1=adjacent cycles only)
spdFiltTag='LowSpd';
%%%%%%%%%%%%%%%%%%%%%%%%

hFig=gra_multiplot(length(SD.selCell),length(SD.selTrial)*2);
ud=get(hFig,'userdata'); hAx=ud.axes_handle_array;
fontspec = {'fontunits','normalized','fontsize', 0.15,'units','normalized'};
axCol=1;

data=evalin('base',SD.selData{1});
for ii=SD.selTrial
    axRow=1;
    % Get best theta channel + theta frequency %
    chInd=data.trials(ii).user.(['bestThetaCh' spdFiltTag]);
    thetaFreq=data.trials(ii).user.(['thetaFreq' spdFiltTag]);
    % Filter the LFP to theta frequency and extract power and phase using Hilbert transform %
    [ freqEEG, phaseEEG, ampEEG ] = eeg_instfrequency(data.trials(ii).eeg(chInd), [thetaFreq-bandWidth thetaFreq+bandWidth]);
    % Calculate power and speed per theta cycle, and exclude bad cycles %
    [cycle, powerPerCycle, speedPerCycle] = eeg_powerspeedpercycle(phaseEEG, ampEEG, data.trials(ii).speed);
    powerThrPctiles= prctile(powerPerCycle,powerThr) ;
    validCycle=find( powerPerCycle>=powerThrPctiles(1) & powerPerCycle<=powerThrPctiles(2) & speedPerCycle>=speedLimit);

    for jj=SD.selCell
        [spikePhase, intrinsicFreq] = spk_intrfreqhilbert(data.trials(ii).cells(jj).st,phaseEEG,cycle,validCycle,maxElapsedCycles,data.trials(ii).eeg(chInd).sample_rate);
        % Histogram of Phase %
        hist(hAx(axRow,axCol*2-1),spikePhase,60);
        set(hAx(axRow,axCol*2-1),'xlim',[0 2*pi]);
        shading(hAx(axRow,axCol*2-1),'flat');
        % Histogram Intrinsic Frequency %
        c=histc(intrinsicFreq,0:0.2:14);       
        bar(hAx(axRow,axCol*2),c,'histc');
        set(hAx(axRow,axCol*2),'xlim',[0 70],'xtick',20:5:50,'xticklabel',4:10);
        shading(hAx(axRow,axCol*2),'flat');
        sumStat=nanmedian(intrinsicFreq);
        text('parent',hAx(axRow,axCol*2),'string',num2str(sumStat,'%3.1f'),'position',[1 1],'verticalalignment','top','horizontalalignment','right',fontspec{1:end});
        line('parent',hAx(axRow,axCol*2),'xdata',repmat(sumStat*5,1,2),'ydata',get(hAx(axRow,axCol*2),'ylim'),'color','g','linewidth',1);
        %%%
        axRow=axRow+1;
    end
    axCol=axCol+1;
end