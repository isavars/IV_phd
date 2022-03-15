function [] = gra_eeg(SD, varargin)
% Plot eeg, with UI for changing scale and section plotted.
%
%       [] = gra_eeg(gss);                To initialise new
%       [] = gra_eeg(gss,'refresh',gcf);  As uicallback

eegInd = 1;

if isempty(varargin)
    %%% First call to function: Initialise figure, axes + UIcontrols %%%
    data = evalin('base', SD.selData{1});
    hFig = gra_multiplot(length(SD.selTrial), 1, 'plotsize', [12 1.5], 'graphborder', [1,1.5,1,1]);
    ud = get(hFig, 'userdata');
    genUISpec = {'units', 'normalized','callback','gra_eeg(gss,''refresh'',gcf);'};
    % Speed filter UI controls (one group for whole figure) %
    axPos = get(ud.axes_handle_array(1),'position');
    speedCheckPos = [0.02 axPos(2)+axPos(4)*5/4 0.08 axPos(4)/5];
    ud.speedCheck = uicontrol(genUISpec{1:end},'style','checkbox','string','Speed Filter','position', speedCheckPos,'tag','speedCheck');
    speedLimitPos = [0.02 axPos(2)+axPos(4)*4/4 0.08 axPos(4)/5];
    ud.speedLimit = uicontrol(genUISpec{1:end},'style','edit','string','2.5','position',speedLimitPos,'tag','speedFilter');
    for ii=1:length(SD.selTrial)
        if isstruct(data.trials(1).eeg)
            % New SCAn data %
            for jj=1:length(data.trials(SD.selTrial(ii)).eeg)
                ud.eeg_data{ii,jj} = (double(data.trials(SD.selTrial(ii)).eeg(jj).eeg) ./ 128) .* data.trials(SD.selTrial(ii)).eeg(jj).scalemax;
                ud.eeg_sample_rate(ii,jj) = data.trials(SD.selTrial(ii)).eeg(jj).sample_rate;
                if length(SD.selCell)==1
                    ud.st{ii}=ceil( data.trials(SD.selTrial(ii)).cells(SD.selCell).st .* data.trials(SD.selTrial(ii)).eeg(jj).sample_rate);
                else
                    ud.st{ii}=[];
                end
                ud.axlim{ii,jj} = [-data.trials(SD.selTrial(ii)).eeg(eegInd).scalemax data.trials(SD.selTrial(ii)).eeg(jj).scalemax];
                ud.selEEGStr{ii}{jj} = ['EEG ' num2str(jj)];
            end
            ud.nEEG(ii)= jj;
            ud.speed{ii}=data.trials( SD.selTrial(ii) ).speed;
            ud.pos_sample_rate(ii)=data.trials( SD.selTrial(ii) ).sample_rate;
            ud.speedFilter{ii}=1:length(ud.eeg_data{ii,jj});
        else
            % Old SCAn data
            ud.eeg_data{ii,1} = data.trials(SD.selTrial(ii)).eeg - 128; % Old scan always from dat files: 0-256
            ud.eeg_sample_rate(ii,1) = data.trials(SD.selTrial(ii)).eeg_sample_rate;
            ud.axlim{ii,1} = [-127 128];
            ud.selEEGStr{ii}{1} = 'EEG 1';
            ud.nEEG(ii)=1;
        end
        ud.hPlotEEG(ii)=nan;    % These will store line handles
        ud.hPlotFilt(ii)=nan;    %
        %%% UI controls %%%
        axPos = get(ud.axes_handle_array(ii),'position');
        xPos1stColUI = 0.85;
        xPos2ndColUI = 0.925;
        % Time limits string box %
        strPos = [xPos1stColUI axPos(2)+(axPos(4)/3) 0.06 axPos(4)/3];
        ud.time_box(ii) = uicontrol(genUISpec{1:end},'style', 'edit', 'string', '0-3', 'position', strPos,'tag','axes');
        % Forward and backward arrows %
        bPos = [xPos1stColUI axPos(2)+(axPos(4)/3)*2 0.03 axPos(4)/3];
        ud.back_button = uicontrol(genUISpec{1:end},'style','push','string','<','position', bPos,'tag',['moveBack' num2str(ii)]);
        fPos = [xPos1stColUI+0.03 axPos(2)+(axPos(4)/3)*2 0.03 axPos(4)/3];
        ud.frwd_button = uicontrol(genUISpec{1:end},'style','push','string','>','position', fPos,'tag',['moveFrwd' num2str(ii)]);
        % EEG Selection Box %
        selPos = [xPos1stColUI axPos(2) 0.06 axPos(4)/4];
        ud.sel_menu(ii) = uicontrol(genUISpec{1:end},'style','popupmenu','string',ud.selEEGStr{ii},'position', selPos,'value',1,'tag','selMenu');
        % Filtering: on/off checkbox %
        filtCheckPos = [xPos2ndColUI axPos(2)+(axPos(4)/4)*3 0.06 axPos(4)/4];
        ud.filtCheck(ii) = uicontrol(genUISpec{1:end},'style','checkbox','string','filter','position', filtCheckPos,'tag',['filtCheck' num2str(ii)]);
        % Filtering: band centre %
        filtFreqPos = [xPos2ndColUI axPos(2)+(axPos(4)/3) 0.025 axPos(4)/4];
        ud.filtFreq(ii) = uicontrol(genUISpec{1:end},'style','edit','string','8','position', filtFreqPos,'tag',['filtFreq' num2str(ii)]);
        % Filtering: band width %
        filtBandPos = [xPos2ndColUI+0.04 axPos(2)+(axPos(4)/3) 0.015 axPos(4)/4];
        ud.filtBand(ii) = uicontrol(genUISpec{1:end},'style','edit','string','3','position', filtBandPos,'tag',['filtBand' num2str(ii)]);
        filtBandLabelPos = [xPos2ndColUI+0.03 axPos(2)+(axPos(4)/3)+(axPos(4)/10) 0.005 axPos(4)/6];
        uicontrol(genUISpec{1:end},'style','text','string',char(177),'position', filtBandLabelPos);
        % Filtering: power per cycle limits %
        powLimLowPos = [xPos2ndColUI axPos(2) 0.025 axPos(4)/4];
        ud.powLimLow(ii) = uicontrol(genUISpec{1:end},'style','edit','string','','position', powLimLowPos,'tag',['powLimLow' num2str(ii)]);
        powLimHighPos = [xPos2ndColUI+0.03 axPos(2) 0.025 axPos(4)/4];
        ud.powLimHigh(ii) = uicontrol(genUISpec{1:end},'style','edit','string','','position', powLimHighPos,'tag',['powLimHigh' num2str(ii)]);

    end
    ud.hFig = hFig;
    ud.trialNames = getdatanames(data,'trial','long',SD.selTrial);
    gra_multilabel(ud.hFig,'row',ud.trialNames);
    set(hFig,'userdata',ud);
else
    %%% Function called as figure callback %%%
    ud = get(varargin{2}, 'userdata');
end


%%% Draw into figure %%%
for ii=1:length(ud.axes_handle_array)
    eegInd = get(ud.sel_menu(ii), 'value');
    tag = get(gcbo,'tag');
    % Have the move backward or forward buttons been pushed? %
    if ~isempty(tag) && strcmp(tag(1:4),'move')
        if str2double(tag(9:end))==ii
            % If have been pushed, change the values in string box %
            [s1 s2] = strtok(get(ud.time_box(ii), 'string'), double('-'));   s2=s2(2:end);
            if strcmp(tag(5:8),'Back')
                newTime = [str2double(s1) str2double(s2)] - (str2double(s2)-str2double(s1));
            else
                newTime = [str2double(s1) str2double(s2)] + (str2double(s2)-str2double(s1));
            end
            set(ud.time_box(ii),'string',[num2str(newTime(1)) '-' num2str(newTime(2))]);
        end
    end
    % Convert time in string boxes to indexes into EEG %
    [ind1sec ind2sec] = strtok(get(ud.time_box(ii), 'string'), double('-'));
    ind1sec = str2double(ind1sec);   ind2sec = str2double(ind2sec(2:end));
    ind1 = round(ind1sec*ud.eeg_sample_rate(ii,eegInd));   ind2 = round(ind2sec*ud.eeg_sample_rate(ii,eegInd));
    if ind1==0;   ind1 = 1;   end;
    % Re-Calculate speed filter, if necessary %
    if ~isempty(tag) && strcmp(tag(1:4),'spee')
        if get(ud.speedCheck,'value')
            ud.speedFilter{ii} = eeg_posind2eegind(  pos_speedfilter(ud.speed{ii}, 'sample_rate',ud.pos_sample_rate(ii), 'min',str2double(get(ud.speedLimit,'string')) )  );   
        else
            ud.speedFilter{ii}=1:length(ud.eeg_data{ii,eegInd});
        end
    end
    % Plot %
    axes(ud.axes_handle_array(ii));
    if ~isnan(ud.hPlotEEG(ii));  delete(ud.hPlotEEG(ii));  ud.hPlotEEG(ii)=NaN;  end  %
    if ~isnan(ud.hPlotFilt(ii));  delete(ud.hPlotFilt(ii));  ud.hPlotFilt(ii)=NaN;  end  % Delete existing plots
    hold on;
    xData=ind1:ind2;
    xData( ~ismember(xData, ud.speedFilter{ii}) ) = NaN;
    ud.hPlotEEG(ii) = plot(xData, ud.eeg_data{ii,eegInd}(ind1:ind2), 'b-');        % Plot eeg
    plot([ind1 ind2], [0 0], 'k:');                            % Plot 0 line
    % Second marker lines % 
    xtickvect = round((ind1sec:ind2sec).*ud.eeg_sample_rate(ii,eegInd));
    for jj=1:length(xtickvect)
        hold on; plot([xtickvect(jj) xtickvect(jj)], [0 256], 'k:');
    end
    set(gca,'xlim',[ind1 ind2],'ylim',ud.axlim{ii,eegInd},'ytick',[ud.axlim{ii,eegInd}(1) 0 ud.axlim{ii,eegInd}(2)],'yticklabel',...
        {num2str(ud.axlim{ii,eegInd}(1),'%3.0f') '0' num2str(ud.axlim{ii,eegInd}(2),'%3.0f')},'xtick',xtickvect,'xticklabel',ind1sec:ind2sec);
    % Cell spikes %
    if ~isempty(ud.st{ii})
        plotSt=ud.st{ii}(ud.st{ii}>ind1 & ud.st{ii}<ind2);
        for jj=1:length(plotSt)
            plot([plotSt(jj) plotSt(jj)],ud.axlim{ii,eegInd},'k-');
        end
    end
    % Filtered EEG signal %
    if get(ud.filtCheck(ii),'value')
        sampleRate = ud.eeg_sample_rate(ii,eegInd);
        lowPass=str2double(get(ud.filtFreq(ii),'string')) - str2double(get(ud.filtBand(ii),'string'));
        highPass=str2double(get(ud.filtFreq(ii),'string')) + str2double(get(ud.filtBand(ii),'string'));
        window = blackman(round(sampleRate)+1);
        EEGFilter = fir1(round(sampleRate), [lowPass highPass]./(sampleRate/2) ,window);
        % Get data to filter (that shown in plot, plus one window length backwards and forwards) %
        filtDataInd = (ind1-length(window)) : (ind2+length(window));    % Index into whole EEG data, gives data to filter
        filtPlotInd = [false(1,length(window)) true(1,length(ind1:ind2)) false(1,length(window))];  % Logical index into filtered data, gives data in the plot window
        overHangInd = filtDataInd<1 | filtDataInd>length(ud.eeg_data{ii,eegInd});   % Check if filtering window overhangs real data ..
        filtDataInd = filtDataInd(~overHangInd);                                    %  .. if so, remove overhang
        filtPlotInd = filtPlotInd(~overHangInd);                                    %  ..
        filtEEG = filtfilt(EEGFilter, 1, ud.eeg_data{ii,eegInd}(filtDataInd));
        if ~isempty(get(ud.powLimLow(ii),'string')) || ~isempty(get(ud.powLimHigh(ii),'string'))
            filtEEG = removeLowPowerCycles(filtEEG, str2double(get(ud.powLimLow(ii),'string')) , str2double(get(ud.powLimHigh(ii),'string')) );
        end
        ud.hPlotFilt(ii) = plot(ind1:ind2,filtEEG(filtPlotInd),'r-');
    end
    hold off;
end
set(ud.hFig,'userdata',ud);
    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rtn]=removeLowPowerCycles(filtEEG,limLow,limHigh)
% Input = filtered EEG signal
% Output = filtered EEG signal with low power cycles set to NaN
analyticEEG = hilbert(filtEEG);
instPhase=angle(analyticEEG);   %Phase (radians) is the angle of the complex variable at each time point
instPhase=unwrap(instPhase);    %Adds 2pi if necessary to render plot of phase smooth
instPower=(abs(analyticEEG)).^2;       %Modulus of analytic function i.e. instantaneous amplitude. Square for amplitude > power
% Split the signal into discrete cycles %
phaseDiff = diff(mod(instPhase, 2*pi)); %Change in phase between each pt - should be a saw tooth
phaseJump = find(phaseDiff<0); %Find point at which phase goes from 2pi to 0
% Remove low power cycles %
nCycles=length(phaseJump)-1;
rtn=filtEEG;
if isnan(limLow);  limLow=0;  end                   % NaNs arise when text inputs are empty
if isnan(limHigh);  limHigh=max(instPower); end     %  ..
for nn = 1:nCycles
    ind=phaseJump(nn):phaseJump(nn+1);
    powerPerCycle=mean(instPower(ind));
    if powerPerCycle<limLow || powerPerCycle>limHigh
        rtn(ind)=nan;
    end
end
    
    
    
    
    
    