
function [h_fig] = gra_allplots(varargin)

% Plot the following for every cell in a given trial: Place map, Dir Map,
% Waveform, Autocorrelogram (20ms) and Autocorrelogram (500ms).
%
%       hFig=gra_allplots(SD)
%           or
%       hFig=gra_allplots(data,cellIndex,trialIndex)
%
% In GUI-called mode, plots cells only for one trial (first selected).
% In command-line mode, if length(trialIndex)=1, plots all cells for that trial. 
% Else, if length(trialIndex)=length(cellIndex), will plot cell(ii) in
% trial(ii) etc.

if length(varargin)==1
    SD = varargin{1};
    data = evalin('base', SD.selData{1});
    cellInd = SD.selCell;
    trInd = SD.selTrial;
    mapInd = SD.selMap;
else
    data = varargin{1};
    cellInd = varargin{2};
    trInd = varargin{3};
    mapInd = 1;
end

     
% mapInd = rates_mapmatch(data.rate_maps, 'bin', 8, 'smooth', 5, 'mode', 'rate');
dirMaps = rates_main(data, rates_params('space','dir','bin',6));   % Make dir maps here, rather than demand them made alread
% Construct unsmoothed maps %
prmsTemp=data.rate_maps(mapInd).params;
prmsDef=rates_params('default');                                    %
f=fieldnames(prmsTemp);                                             % These lines are needed to protect against old maps, which do not have all the params fields that RATES_MAIN is expecting.
for ii=1:length(f);   prmsDef.(f{ii}) = prmsTemp.(f{ii});  end      %
prmsDef.smooth=1;
unSmoothedMaps = rates_main(data,prmsDef);

% Set up figure %
h_fig = gra_multiplot(length(cellInd), 7, 'graphborder', repmat(0.7, [1 4]), 'figborder', [1 0.5 0.5 0.5]);
set(h_fig, 'defaulttextfontname', 'arial', 'numbertitle', 'off');
set(h_fig,'inverthardcopy','off','color','white');      % Workaround for Matlab bug that prints white as text black with 'axis off'.
% Menus %
set(h_fig, 'menubar', 'none');   % Turn off matlab default.
uimenu('label', 'Save as eps ..', 'callback', 'scan_base(''saveEps'',gcf);');
uimenu('label', 'Print', 'callback', 'printdlg(gcf);');

% Draw Plots %
temp = get(h_fig, 'userdata');   hArray = temp.axes_handle_array;
axCountR=1;
for ii=cellInd
    %%%% (1) Place Map %
    axes(hArray(axCountR,1));
    gra_plotmap(data.rate_maps(mapInd).maps{trInd, ii});
    
    %%%% (2) Spatial autocorr + grid props %
    % Gridness method is 'Barry 2007'
    % 'The six peaks surrounding the central peak on the autocorrelogram were considered to be the local maxima closest to,
    % but excluding, the central peak. The extent of each peak was defined as the contiguous set of bins around the peak with 
    % a value greater than half the value of the peak bin.'
    axes(hArray(axCountR,2));
    map=unSmoothedMaps{trInd,ii};
    ac = map_crosscorr(map,map,'smooth', 10, 2.25);
    [gr, grProps] = map_gridprops(ac,'peakMode','point','corrThr',0,'radius','hh');
    ac(isnan(ac)) = min(min(ac)) - ( (max(max(ac))-min(min(ac)))*0.1 ); % Auto-scale unvisited bins to first value.
    colormap(gra_tintcolormap(10));
    imagesc(ac);  axis off;  axis equal;   
    % Mark close peaks %
    line('xdata', grProps.closestPeaksCoord(:,1), 'ydata', grProps.closestPeaksCoord(:,2), 'linestyle', 'none', 'marker', 's', 'color', 'k');
    % Mark gridness calculation boundary %
    [x,y] = pol2cart(0:pi/15:2*pi, repmat(grProps.maxDistFromCentre,1,31));
    line('xdata', x+(size(ac,2)/2), 'ydata', y+(size(ac,1)/2), 'marker', 'none', 'color', [0.0 0.0 0.0]);
    % Gridness text %
    fontspec = {'fontunits','normalized','fontsize',0.15,'color','k','interpreter','tex'};
    text('string', num2str(gr, '%3.2f'),'position', [0 1],'units','normalized','verticalalignment','cap','horizontalalignment','left', fontspec{1:end});
    
    %%%% (3) Dir map %
    axes(hArray(axCountR,3));
    gra_plotmap(dirMaps{trInd, ii});
    %%%% (4)  Waveform %
    axes(hArray(axCountR,4));
    gra_plotwf(data.trials(trInd).cells(ii).wf_means,data.trials(trInd).cells(ii).scalemax,'max',gca);
    %%%% (5)  Temporal AC, short range %
    axes(hArray(axCountR,5));
    spk_crosscorr(data.trials(trInd).cells(ii).st,data.trials(trInd).cells(ii).st,0.5,20,data.trials(trInd).dur,'plot',gca);
    %%%% (6)  Temporal AC, long range %
    axes(hArray(axCountR,6));
    spk_crosscorr(data.trials(trInd).cells(ii).st,data.trials(trInd).cells(ii).st,10,500,data.trials(trInd).dur,'plot',gca);
    
    %%%% (7) All 4 waveforms %%%%
    axPos=get(hArray(axCountR,7),'position');
    delete(hArray(axCountR,7));
    hNewAx(1) = axes('position',[axPos(1:2) axPos(3:4)./2]);
    hNewAx(2) = axes('position',[((axPos(3)/2) + axPos(1)) axPos(2) axPos(3:4)./2]);
    hNewAx(3) = axes('position',[axPos(1) ((axPos(4)/2) + axPos(2))  axPos(3:4)./2]);
    hNewAx(4) = axes('position',[((axPos(3)/2) + axPos(1)) ((axPos(4)/2) + axPos(2))  axPos(3:4)./2]);
    H=gra_plotwf(data.trials(trInd).cells(ii).wf_means,data.trials(trInd).cells(ii).scalemax,'all',hNewAx);
    for jj=[2 4]
        posTemp = get(H.handles(jj).yAxis_text,'position');
        posTemp(1) = 63;
        set(H.handles(jj).yAxis_text,'position',posTemp,'horizontalalignment','right');
    end
    for jj=1:4
        set(H.handles(jj).yAxis_text,'fontsize',0.2);
    end
    delete(H.handles(H.maxChannel).props_text);
    
    % Finished %
    axCountR=axCountR+1;
end




% Label (Only if access to data name - GUI calling mode) %
if length(varargin)==1
    cStr=getdatanames(data, 'cell', 'short', cellInd);   trStr=getdatanames(data, 'trial', 'long', trInd);
    for ii=1:length(trInd);  lStr{ii}(1)=cStr(ii);  lStr{ii}(2)=trStr(ii);   end
    hRowLab = gra_multilabel(h_fig, 'row', lStr);
    set(hRowLab,'fontsize',0.4);
    hTitle = gra_multilabel(h_fig, 'title', SD.selData{1});
    set(h_fig,'name', SD.selData{1});
end



