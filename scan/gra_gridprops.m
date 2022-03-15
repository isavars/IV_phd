function [hFig] = gra_gridprops(input1,gridScoreType,varargin)
% Draw rate map auto correlograms, and label with gridness properties.
% Can call from SCAn GUI, or from command line:
% 
%         gra_gridprops(SD,'gridScoreType');    (As called form SCAn GUI)
%               or
%         hFig = gra_gridprops(data,'gridScoreType',trialIndex,cellIndex,rateMapInd);      (To call from code, with data struct)

% This is not in use ..
%               or
%         hFig = gra_gridprops([],'gridScoreType',trialIndex,cellIndex,rateMaps,ACMaps,binsize,trialLabels,cellLabels);    (To call from code, without data struct)

if isstruct(input1) && isfield(input1,'selData')
    % If called by SCAn GUI %
    data = evalin('base', input1.selData{1});
    trInd=input1.selTrial;
    cellInd=input1.selCell;
    rateMaps = data.rate_maps( input1.selMap ).maps;
    rateMapInd = input1.selMap;
    binsize = (1/data.trials(1).ppm) * data.rate_maps( rateMapInd ).params.bin * 100;
    mapsUnsmoothed = length(data.rate_maps(rateMapInd).params.smooth)==1  && data.rate_maps(rateMapInd).params.smooth==1;
    
elseif isstruct(input1) && isfield(input1,'trials')
    % If called in code, but using SCAn data struct %
    data = input1;
    trInd=varargin{1};
    cellInd=varargin{2};
    rateMapInd = varargin{3};
    rateMaps = data.rate_maps(rateMapInd).maps;
    binsize = (1/data.trials(1).ppm) * data.rate_maps( rateMapInd ).params.bin * 100;
    mapsUnsmoothed = length(data.rate_maps(rateMapInd).params.smooth)==1  && data.rate_maps(rateMapInd).params.smooth==1;
    
    % NOT IN USE
% elseif isempty(input1)
%     % If called in code, but without SCAn data struct %
%     trInd=varargin{1};
%     cellInd=varargin{2};
%     rateMaps = varargin{3};
%     ACMaps = varargin{4};
%     binsize  = varargin{5};
%     mapsUnsmoothed = 0;     % In this case, we can't know if the maps are smoothed or not - will have to trust the user to supply correct map type.
    
end

% Check if we need to construct unsmoothed maps (for Barry 2007 gridness) (do this by getting maps params struct, setting smooth=1, and remaking maps) %
if ~mapsUnsmoothed && any(strcmp(gridScoreType,{'Barry 2007','Lever 2012'})); 
    disp('Constructing unsmoothed rate maps .. ');
    prmsTemp=data.rate_maps(rateMapInd).params;
    prmsDef=rates_params('default');                                    %
    f=fieldnames(prmsTemp);                                             % These lines are needed to protect against old maps, which do not have all the params fields that RATES_MAIN is expecting.
    for ii=1:length(f);   prmsDef.(f{ii}) = prmsTemp.(f{ii});  end      %
    prmsDef.smooth=1;
    ACMaps=rates_main(data,prmsDef);
else
    ACMaps=rateMaps;
end

hFig = gra_multiplot(length(cellInd), length(trInd)*2, 'graphborder', repmat(0.2, [1 4]));    ud = get(hFig, 'userdata');
trial_count = 1;
gridPropsVerbose = 0;
for ii=trInd 
    cell_count = 1;
    for jj=cellInd
        map = ACMaps{ii,jj};
        if nanmax(nanmax(map))==0
            cell_count=cell_count+1;
            continue
        end

        % Calculate gridness %
        switch gridScoreType
            case 'Sargolini 2006'
                % peaks defined as 100 or more contiguous pixels of 1.5 × 1.5 cm2 exceeding an arbitrary threshold, in most cases r=0.10
                ac = map_crosscorr(map,map);
                minPeakArea = 225/(binsize^2);
                [gr, grProps] = map_gridprops(ac,'peakMode','area','areaThr',minPeakArea,'corrThr',0.1,'radius','fieldExtent','verbose',gridPropsVerbose);
            case 'Barry 2007'
                % 'The six peaks surrounding the central peak on the autocorrelogram were considered to be the local maxima closest to,
                % but excluding, the central peak. The extent of each peak was defined as the contiguous set of bins around the peak with 
                % a value greater than half the value of the peak bin.'
                ac = map_crosscorr(map,map,'smooth');
                [gr, grProps] = map_gridprops(ac,'peakMode','point','corrThr',[],'radius','fieldExtent','verbose',gridPropsVerbose);
            case 'Wills 2010'
                ac = map_crosscorr(map,map,'smooth');
                [gr, grProps] = map_gridprops(ac,'peakMode','point','corrThr',0,'radius','est','verbose',gridPropsVerbose);
            case 'Watershed'
                ac = map_crosscorr(map,map,'smooth');
                [gr, grProps] = map_gridprops(ac,'peakMode','point','corrThr',0,'radius','fieldExtent','verbose',gridPropsVerbose,'fieldExtentMethod','watershed');
            case 'Lever 2012'
                % Similar to Barry 2007 method, peaks must be r>0, larger
                % smoothing kernel for AC
                ac = map_crosscorr(map,map,'smooth', 10, 2.25);
                [gr, grProps] = map_gridprops(ac,'peakMode','point','corrThr',0,'radius','fieldExtent');
                
%             case 'Langston 2010'
%                  [gr]=dev_langstongridness(ac,binsize,10:5:40,'trondheim');
        end
        gridPropsVerbose = 0;        
            
        % Plot AC and print grid properties %
        axes(ud.axes_handle_array(cell_count,trial_count+1));
        ac(isnan(ac)) = min(min(ac)) - ( (max(max(ac))-min(min(ac)))*0.1 ); % Auto-scale unvisited bins to first value.
        colormap(gra_tintcolormap(10));
        imagesc(ac);   
        axis off;  axis square;   axis equal;
        if ~isnan(grProps.closestPeaksCoord)
            % Mark close peaks %
            line('xdata', grProps.closestPeaksCoord(:,1), 'ydata', grProps.closestPeaksCoord(:,2), 'linestyle', 'none', 'marker', 's', 'color', 'k');
            % Mark gridness calculation boundary %
            [x,y] = pol2cart(0:pi/15:2*pi, repmat(grProps.maxDistFromCentre,1,31));
            line('xdata', x+(size(ac,2)/2), 'ydata', y+(size(ac,1)/2), 'marker', 'none', 'color', [0.0 0.0 0.0]);
        end
        % Gridprops text %
        fontspec = {'fontunits','normalized','fontsize',0.1,'color','k','interpreter','tex'};
        text('string', num2str(gr, '%3.2f'),'position', [0 1],'units','normalized','verticalalignment','top','horizontalalignment','left', fontspec{1:end});
        text('string', ['\theta=' num2str(grProps.orientation, '%3.0f')],'units','normalized','position',[1 1],'units','normalized','verticalalignment','top','horizontalalignment', 'right', fontspec{1:end});
        text('string', ['\lambda=' num2str(grProps.waveLength*binsize,'%3.1f')], 'units', 'normalized', 'position', [0 0],'verticalalignment','bottom','horizontalalignment','left',fontspec{1:end});
        
        % Plot Rate Map %
        axes(ud.axes_handle_array(cell_count,trial_count));
        gra_plotmap(rateMaps{ii,jj}); 
        
        cell_count = cell_count + 1;
    end
    trial_count = trial_count + 2;
end
% Label (only if called GUI, otherwise need to self-label %
if ~isempty(data) % In case data geiven as a dummy, only maps passed.
    t = getdatanames(data,'trial','short',trInd);
    tLabel(2:2:length(t)*2) = t;
    gra_multilabel(hFig, 'col', tLabel); 
    gra_multilabel(hFig, 'row', getdatanames(data,'cell','short',cellInd));
    if isfield(input1,'selData')
        gra_multilabel(hFig, 'title', input1.selData{1});
    end
end






