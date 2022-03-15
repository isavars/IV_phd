function [] = spikesOnDirGUI(varargin)
% Plot Head Direction, Angular Velocity, with UI for changing scale and section plotted.
%
%       [] = spikesOnDirGUI();                    To initialise new
%       [] = spikesOnDirGUI(gss,'refresh',gcf);   As uicallback

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET UP GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin)
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get data necessary for UI: dir, ang velocity, spike times %
    SD=gss;
    data = evalin('base', SD.selData{1});
    ud.pos_sample_rate = data.trials(SD.selTrial(1)).sample_rate;
    d = double(data.trials(SD.selTrial(1)).dir);
    ud.AV_data = abs(diff( unwrap( deg2rad(  d  ) ) )) .* ud.pos_sample_rate .* (180/pi);
    ud.AV_data(end+1) = nan;
    % Pre-process the HD data, so as to remove excessive 'wraps' or jumps around 0/360 %
    dDiff = diff(d);
    jumpList = find(  abs(dDiff) > 300  );  % List of jumps across 0/360
    jumpPol = ones(size(jumpList));              % Get the polarity of each jump
    jumpPol(  dDiff(jumpList) < -300  ) = -1;    %  ..
    ii=1;
    while ii<=(length(jumpList)-1)
        if jumpPol(ii)~=jumpPol(ii+1)  &&   ((jumpList(ii+1)-jumpList(ii))/ud.pos_sample_rate) < 10 % If we are looking at jump that goes forwards and then backwards, within 10 s ..
            jumpInd = (jumpList(ii)+1):jumpList(ii+1);
            dWithinJump = d(  jumpInd  );
            if jumpPol(ii)==-1 && nanmax(dWithinJump)<=90
                d( jumpInd ) = d(jumpInd) + 360;   ii=ii+1;
            elseif jumpPol(ii)==1 && nanmin(dWithinJump)>=270
                d( jumpInd ) = d(jumpInd) - 360;   ii=ii+1;
            end
        end
        ii=ii+1;
    end     
    d(  abs(diff(d)) > 300  ) = NaN;   % Set the pre-jump bin of remaining jumps to NaN, so that the 'jumping line' isn't plotted.
    ud.HD_data = d;
    % Spike times: store them all here, then up to 4 can be selected for plotting via GUI (see below) %
    for ii=1:length(data.trials(SD.selTrial(1)).cells)
        trData = data.trials(SD.selTrial(1));
        ud.ST_all{ii}=ceil( trData.cells(ii).st .* trData.sample_rate);
        ud.cellNamesAll{ii+1} = ['T' num2str(trData.cells(ii).tet) ' C' num2str(trData.cells(ii).cellnum)];
    end
    ud.cellNamesAll{1} = 'Not Active';  % First option in list is not to show any cell data, for that slot.
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create figure %
    % We want a bigger axis for the HD signal than everything else. As a quick fix, call GRA_MULTIPLOT with more 
    % rows than we need, then fuse 3 of these plots into one for the HD data axis.
    hFig = gra_multiplot(8, 1, 'plotsize', [12 1.5], 'graphborder', [1,1.5,1,1]);   axArray = getappdata(hFig,'axesHandles'); 
    hTemp = axArray(2);   
    axTopPos = hTemp.Position(2) + hTemp.Position(4);
    axPosNew = get(axArray(4),'position');
    axPosNew(4) = axTopPos - axPosNew(2);
    set(axArray(4),'position',axPosNew);
    delete(axArray(2:3));
    
    % Assign Axes handles for graphs to GUI struct %
    ud.hAxisAV=axArray(1);
    ud.hAxisHD=axArray(4);
    for ii=1:4
        ud.hAxisSpkHist(ii) = axArray(ii+4);
    end
    % Pre-assign storage for handles for graphics objects %
    ud.hPlotHD=nan;          % These will store line handles
    ud.hPlotAV=nan;          %
    ud.hPlotSpikes=nan(1,4);   % These will store the handles for the point markers (also actually line objects) which show spikes on HD.
    ud.hPlotSpkHist=nan(1,4);  % These are for the rate histgram objects
    
    %%% UI controls %%%
    genUISpec = {'units', 'normalized','callback','spikesOnDirGUI(gss,''refresh'',gcf);'};
    axPos = get(ud.hAxisAV,'position');
    xPos1stColUI = 0.85;
%     xPos2ndColUI = 0.925;
    % Time limits string box %
    strPos = [xPos1stColUI axPos(2)+(axPos(4)/3) 0.06 axPos(4)/3];
    ud.time_box = uicontrol(genUISpec{1:end},'style', 'edit', 'string', '0-60', 'position', strPos,'tag','');
    % Forward and backward arrows %
    bPos = [xPos1stColUI axPos(2)+(axPos(4)/3)*2 0.03 axPos(4)/3];
    ud.back_button = uicontrol(genUISpec{1:end},'style','push','string','<','position', bPos,'tag','moveBack' );
    fPos = [xPos1stColUI+0.03 axPos(2)+(axPos(4)/3)*2 0.03 axPos(4)/3];
    ud.frwd_button = uicontrol(genUISpec{1:end},'style','push','string','>','position', fPos,'tag','moveFrwd');
    % Time bin for rate histogram %
    axPos = get(ud.hAxisSpkHist(1),'position');
    hTimePos = [xPos1stColUI axPos(2)+(axPos(4)/5)*3 0.06 axPos(4)/5];
    ud.spkHistTimeBin = 1.0;  % This entry sets the default time bin for the rate histograms.
    uicontrol('units', 'normalized','style','text','string','Time Div (sec):','position',hTimePos+[0 axPos(4)/5 0 0]);
    ud.hSpkHistTimeBin_strBox = uicontrol(genUISpec{1:end},'style', 'edit', 'string', num2str(ud.spkHistTimeBin), 'position', hTimePos,'tag','hSpkHistTimeBin_strBox');
      
    % Cell Selection Drop-downs %
    for ii=1:4
        vPos = ud.hAxisSpkHist(ii).Position(2);
        selPos = [xPos1stColUI, vPos, 0.06, axPos(4)/4];
        ud.cellDropDown(ii) = uicontrol(genUISpec{1:end},'style','popupmenu','string',ud.cellNamesAll,'position', selPos,'value',1,'tag',['cellDropDown' num2str(ii)]);
    end
    
    % Set the cell selection (first two cells selcted in scan) %
    ud.cellsPlotted = zeros(1,4);  % This field defines which cell data populates the spike hist plots (and also coloured dots on HD line). 0=plot currently empty, plot no data here.
    for ii=1:min([4, length(SD.selCell)]);
        ud.cellsPlotted(ii) = SD.selCell(ii);
        set(ud.cellDropDown(ii),'value',SD.selCell(ii)+1);
    end
     
    % Save the GUI struct to the figure userdata %
    ud.hFig = hFig;
    ud.trialName = getdatanames(data,'trial','long',SD.selTrial(1));
    set(hFig,'userdata',ud);

else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function called as figure callback %
    ud = get(varargin{3}, 'userdata');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POPULATE FIGURE WITH DATA, REFRESH WHEN BUTTON IS PRESSED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dectect GUI interaction %%%
tag = get(gcbo,'tag');
% (1) Have the move backward or forward in time buttons been pushed? 
if ~isempty(tag) && strncmp(tag,'move',4)
    % If have been pushed, change the values in string box - NOTE that the string box is the ultimate determiner of what the time span of the plots is %
    [s1, s2] = strtok(get(ud.time_box, 'string'), double('-'));   s2=s2(2:end);
    if strcmp(tag(5:8),'Back')
        newTime = [str2double(s1) str2double(s2)] - (str2double(s2)-str2double(s1));
    else
        newTime = [str2double(s1) str2double(s2)] + (str2double(s2)-str2double(s1));
    end
    set(ud.time_box,'string',[num2str(newTime(1)) '-' num2str(newTime(2))]);
end
% (2) Have the cell selections been changed  %
if ~isempty(tag) && strncmp(tag,'cellDropDown',12)
    ud.cellsPlotted( str2double(tag(end)) ) = get(gcbo,'value') - 1;  % Cell index is value-1, as first entry in cell list is 'NOT ACTIVE'.
end
% (3) Has the rate histogram time bin changed %
if ~isempty(tag) && strcmp(tag,'hSpkHistTimeBin_strBox')
    ud.spkHistTimeBin = str2double(  get(gcbo,'string')  );
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert time in string boxes to indexes into EEG %
[ind1sec, ind2sec] = strtok(get(ud.time_box, 'string'), double('-'));
ind1sec = str2double(ind1sec);   
ind2sec = str2double(ind2sec(2:end));
ind1 = round(ind1sec*ud.pos_sample_rate);   ind2 = round(ind2sec*ud.pos_sample_rate);
if ind1==0;   ind1 = 1;   end;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the existing plots %
if ishandle(ud.hPlotHD);  delete(ud.hPlotHD);    end  %
if ishandle(ud.hPlotAV);  delete(ud.hPlotAV);    end
for ii=1:4
    if ishandle(ud.hPlotSpikes(ii));   delete(ud.hPlotSpikes(ii));    end
    if ishandle(ud.hPlotSpkHist(ii));   delete(ud.hPlotSpkHist(ii));    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Main Data %
hold(ud.hAxisHD,'on');   hold(ud.hAxisAV,'on');
xData=ind1:ind2;
%     xData( ~ismember(xData, ud.speedFilter{ii}) ) = NaN;
ud.hPlotHD = plot(ud.hAxisHD, xData, ud.HD_data(ind1:ind2), 'k-');       % Plot HD
ud.hPlotAV = plot(ud.hAxisAV, xData, ud.AV_data(ind1:ind2), 'k-');       % Plot AV    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Cell spikes %
cList = 'rbmg';
for ii=1:4  % ii=interator for the 4 active cell slots
    if ud.cellsPlotted(ii)~=0
        % (1) Plot the spikes on the HD curve %
        stForPlot = ud.ST_all{ud.cellsPlotted(ii)};
        stForPlot = stForPlot(stForPlot>ind1 & stForPlot<ind2);
        if ~isempty(stForPlot)  
            ud.hPlotSpikes(ii)=plot(ud.hAxisHD,stForPlot,ud.HD_data(stForPlot),[cList(ii),'d'],'MarkerSize',2);
        end
        % (2) Draw the firing rate histograms %
        ud.hPlotSpkHist(ii)=histogram( ud.hAxisSpkHist(ii), stForPlot, ind1:(ud.spkHistTimeBin*ud.pos_sample_rate):ind2 );
        set(ud.hPlotSpkHist(ii),'FaceColor',cList(ii),'EdgeColor','none');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format Axes %                     
% Calculate the x-axis division markers, and their matching labels (in sec) %
possTimeDivs = [1 2 5 10 20 30 60];   % In s
[~,bestTimeDiv] = min(     abs(    ((ind2sec-ind1sec)./possTimeDivs) - 5        )   );
xTickTimeDiv = possTimeDivs(bestTimeDiv); 
xTickLabel = ind1sec:xTickTimeDiv:ind2sec;
xTickVect = xTickLabel .* ud.pos_sample_rate;
plot(ud.hAxisHD, [ind1 ind2], [0 0], 'k-');    plot(ud.hAxisHD, [ind1 ind2], [360 360], 'k-');
% Set the Axis properties %
set(ud.hAxisHD,'xlim',[ind1 ind2],'ylim',[-90 450],'ytick',-90:90:450,'yticklabel',-90:90:450,'xtick',xTickVect,'xticklabel',xTickLabel,'ygrid','on','xgrid','on');
set(ud.hAxisAV,'xlim',[ind1 ind2],'xtick',xTickVect,'xticklabel',xTickLabel,'xgrid','on','yticklabelmode','auto');
for ii=1:4
    set(ud.hAxisSpkHist,'xlim',[ind1 ind2],'xtick',xTickVect,'xticklabel',xTickLabel,'ygrid','on');
end
ylabel(ud.hAxisAV,{'AHV','Deg/S'});
ylabel(ud.hAxisHD,{'HD','Deg'});
for ii=1:4
    ylabel(ud.hAxisSpkHist(ii),{['Cell ' num2str(ii) ' Rate'],'Hz'});
end

% Save the UI struct back to the figure %
set(ud.hFig,'userdata',ud);
    
 

    
    
    
    
    