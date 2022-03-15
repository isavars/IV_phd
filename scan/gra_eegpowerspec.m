function []=gra_eegpowerspec(SD)
% Plot power spectra of EEG and distribution of power per cycle when
% filtered at peak theta frequency.

%%%%%%%%%%%%%%%%%%%%%%%%
bandWidth=3;            % This is actually half-band width (pass band is peak+-bandWidth)
hfCutOff=25;            % High frequency cut-off for power spectra
speedLimit=5;
medianSpeed=5.78;
speedChunk=0.5;
powerThr=[10 99];        % Inclusion cut-offs for power per cycle (%tiles of powerPerCyle for whole trial). Applied after speed filtering.
%%%%%%%%%%%%%%%%%%%%%%%%

data=evalin('base',SD.selData{1});
hFig=gra_multiplot(length(SD.selTrial),4);
ud=get(hFig,'userdata'); hAx=ud.axes_handle_array;

for ii=1:length(SD.selTrial)
    % Calculate speed filter %
    [posFiltInd,dum1,dum2]=pos_speedfilter(data.trials(SD.selTrial(ii)),'min',speedLimit,'chunk',speedChunk,'median',[]);
    eegFiltInd=eeg_posind2eegind(posFiltInd);
    speedAtEEGSampRate=repmat(data.trials(SD.selTrial(ii)).speed', 5, 1); % col = 1 speed samp, row=1:5 reps of speed samp
    speedAtEEGSampRate=reshape(speedAtEEGSampRate,numel(speedAtEEGSampRate),1); % Reshape to column vector
    % For each EEG channel .. %
    for jj=1:2
        if jj>length(data.trials( SD.selTrial(ii) ).eeg)
            continue
        end
        sr=data.trials(SD.selTrial(ii)).eeg(jj).sample_rate;
        % Calculate and plot powerspec %
        [peakFreq dum dum]=eeg_powerspec( data.trials(SD.selTrial(ii)).eeg(jj).eeg( eegFiltInd ), sr, 'freqRes', 0.01, 'hfCutOff', hfCutOff, 'fig',1,'hAxis',hAx(ii,jj*2-1) );
        % Filter EEG to theta frequency %
        window = blackman(round(sr)+1);
        EEGFilter = fir1(round(sr), [peakFreq-bandWidth peakFreq+bandWidth]./(sr/2) ,window);
        filtEEG = filtfilt(EEGFilter, 1, double(data.trials(SD.selTrial(ii)).eeg(jj).eeg));
        % Calculate distribution of power per cycle %
        analyticEEG = hilbert(filtEEG); %Hilbert transform itself - returns a complex var.
        phaseEEG250hz=angle(analyticEEG); %Phase (radians) is the angle of the complex variable at each time point
        eegPhase=unwrap(phaseEEG250hz); %Adds 2pi if necessary to render plot of phase smooth
        eegAmp=abs(analyticEEG); %Modulus of analytic function i.e. instantaneous amplitude
        [cycle, powerPerCycle, speedPerCycle] = findPowerPerCycle(eegPhase,eegAmp,speedAtEEGSampRate); % This can be changed for new function eeg_powerspeedpercycle, when it is formalised in SCAn    
        powerPerCycle=powerPerCycle(speedPerCycle>=speedLimit);
        % Plot distribution of power per cycle %
        hist(hAx(ii,jj*2),powerPerCycle,50);
        fontspec = {'fontunits','normalized','fontsize', 0.075,'units','normalized'};
        pc=round( prctile(powerPerCycle,powerThr) );
        for kk=1:length(pc);  line('xdata',[pc(kk) pc(kk)],'ydata',get(hAx(ii,jj*2),'ylim'),'parent',hAx(ii,jj*2));  end
        text('parent',hAx(ii,jj*2),'string',num2str(pc(1)),'position',[1 1],'verticalalignment','top','horizontalalignment','right',fontspec{1:end});
        text('parent',hAx(ii,jj*2),'string',num2str(pc(2)),'position',[1 0.9],'verticalalignment','top','horizontalalignment','right',fontspec{1:end}); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cycle, powerPerCycle, speedPerCycle]=findPowerPerCycle(eegPhase,eegAmp,speed)
% Find power per cycle + filter for slow speeds included in cycle
phaseDiff = diff(mod(eegPhase, 2*pi)); %Change in phase between each pt - should be a saw tooth
phaseJump = find(phaseDiff<0); %Find point at which phase goes from 2pi to 0

nCycles=length(phaseJump)-1;
powerPerCycle=zeros(nCycles,1); %Pre allocate for speed
speedPerCycle=zeros(nCycles,1);
cycle=zeros(size(eegPhase)); %Vector same size as eegPhase
eegAmp=eegAmp.^2; %Square to turn amplitude to power

for nn = 1:nCycles
    %Get mean power by taking the mean of already squared amplitude
    speedPerCycle(nn) = mean(speed( phaseJump(nn):phaseJump(nn+1) ));
    powerPerCycle(nn)=mean(eegAmp(phaseJump(nn):phaseJump(nn+1)));
    cycle(phaseJump(nn):phaseJump(nn+1))=nn;
end