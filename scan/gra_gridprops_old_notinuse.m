function [hFig] = gra_gridprops_old_notinuse(input1,varargin)
% Draw rate map auto correlograms, and label with gridness properties.
% Can call from SCAn GUI, or from command line:
% 
%         gra_gridprops(SD)
%               or
%         hFig = gra_gridprops(data,trialIndex,cellIndex,rateMaps,ACMaps,binsize);

if isfield(input1,'selData')
    % If called by SCAn GUI %
    data = evalin('base', input1.selData{1});
    trInd=input1.selTrial;
    cellInd=input1.selCell;
    rateMaps = data.rate_maps( rates_mapmatch(data,input1.selMap) ).maps;
    ACMapInd = rates_mapmatch(data,input1.selMap);
    ACMaps = data.rate_maps( ACMapInd ).maps;
    binsize = (1/data.trials(1).ppm) * data.rate_maps( ACMapInd ).params.bin * 100;
    mapsUnsmoothed = length(data.rate_maps(ACMapInd).params.smooth)==1  && data.rate_maps(ACMapInd).params.smooth==1;
else
    % If called in code %
    data = input1;
    trInd=varargin{1};
    cellInd=varargin{2};
    rateMaps = varargin{3};
    ACMaps = varargin{4};
    binsize  = varargin{5};
    mapsUnsmoothed = 0;
end
    
hFig = gra_multiplot(length(cellInd), length(trInd)*2, 'graphborder', repmat(0.2, [1 4]));
ud = get(hFig, 'userdata');

trial_count = 1;
gridPropsVerbose = 1;
for ii=trInd 
    cell_count = 1;
    for jj=cellInd
        map = ACMaps{ii,jj};
        if nanmax(nanmax(map))==0
            cell_count=cell_count+1;
            continue
        end
        % Plot Rate Map %
        axes(ud.axes_handle_array(cell_count,trial_count));
        gra_plotmap(rateMaps{ii,jj});
        % Plot Autocorrelogram %
        if mapsUnsmoothed
            % If input rate map is unsmoothed, smooth AC %
            ac = map_crosscorr(map,map,'smooth');
            [gr, grProps] = map_gridprops(ac,'peakMode','point','corrThr',0,'radius','est','verbose',gridPropsVerbose);
        else
            % .. otherwise, don't smooth AC %
            ac = map_crosscorr(map,map);
            [gr, grProps] = map_gridprops(ac,'peakMode','point','areaThr',20,'corrThr',0.3,'verbose',gridPropsVerbose);
%             [gr]=dev_langstongridness(ac,binsize,10:5:40,'trondheim');
        end        
        axes(ud.axes_handle_array(cell_count,trial_count+1));
        ac(isnan(ac)) = min(min(ac)) - ( (max(max(ac))-min(min(ac)))*0.1 ); % Auto-scale unvisited bins to first value.
%         ac(peakMask)=max(max(ac));  % Mark whole extent of centre peak (for field size)
        colormap(gra_tintcolormap(10));
        imagesc(ac);   axis square;   axis off;
%         if ~isnan(pks)
%             % Mark close peaks %
%             line('xdata', grProps.closestPeaksCoord(:,1), 'ydata', grProps.closestPeaksCoord(:,2), 'linestyle', 'none', 'marker', 's', 'color', 'k');
%             % Mark gridness calculation boundary %
%             [x,y] = pol2cart(0:pi/15:2*pi, repmat(grProps.maxDistFromCentre,1,31));
%             line('xdata', x+(size(ac,2)/2), 'ydata', y+(size(ac,2)/2), 'marker', 'none', 'color', [0.0 0.0 0.0]);
%         end
        % Gridprops text %
        if gr>0;   fontColor = 'k';   else   fontColor = 'g';   end
        fontspec = {'fontunits', 'normalized', 'fontsize', 0.12};
        text('string', num2str(gr, '%3.2f'), 'position', [0 1], 'units', 'normalized', ...
             'verticalalignment', 'top', 'horizontalalignment', 'left', fontspec{1:end},'color',fontColor);
        text('string', num2str(grProps.orientation, '%3.1f'), 'position', [1 1], 'units', 'normalized', ...
             'verticalalignment', 'top', 'horizontalalignment', 'right', fontspec{1:end});
%         text('string', num2str(grProps.orientation, '%3.1f'), 'position', [1 0], 'units', 'normalized', ...
%              'verticalalignment', 'bottom', 'horizontalalignment', 'right', fontspec{1:end});
%         text('string', num2str(wl*binsize,'%3.1f'), 'position', [0 0], 'units', 'normalized', ...
%              'verticalalignment', 'bottom', 'horizontalalignment', 'left', fontspec{1:end});
        cell_count = cell_count + 1;
        gridPropsVerbose = 0;
    end
    trial_count = trial_count + 2;
end
% Label (only if called GUI, otherwise need to self-label %
if ~isempty(data) % Incase data geiven as a dummy, only maps passed.
    t = getdatanames(data,'trial','short',trInd);
    tLabel(2:2:length(t)*2) = t;
    gra_multilabel(hFig, 'col', tLabel); 
    if isfield(input1,'selData')
        gra_multilabel(hFig, 'row', getdatanames(data,'cell','short',cellInd));
        gra_multilabel(hFig, 'title', input1.selData{1});
    end
end






