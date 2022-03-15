function gra_spatialcrosscorr(SD,varargin)
% Spatial cross-correlograms between rate maps.

data = evalin('base', SD.selData{1});
trInd=SD.selTrial;
cellInd=SD.selCell;
mapInd=rates_mapmatch(data,SD.selMap);

prms.trialMode='across';
% prms.acBin=5;
% prms.acWin=100;
prms.addSpatialMap=0;
for ii=1:2:length(varargin)
    prms.(varargin{ii}) = varargin{ii+1};
end

switch prms.trialMode
    case 'across'
        hFig = gra_multiplot(length(cellInd),3,'graphborder',repmat(0.2, [1 4]));   
        ud=get(hFig, 'userdata');   ax=ud.axes_handle_array;
        for ii=1:length(cellInd)
            mapA=data.rate_maps(mapInd).maps{trInd(1),cellInd(ii)};
            mapB=data.rate_maps(mapInd).maps{trInd(2),cellInd(ii)};
            plotCorrelogram(mapA,mapA,ax(ii,1));    % Autocorr for the two trials in the first two columns ..
            plotCorrelogram(mapB,mapB,ax(ii,2));    % 
            plotCorrelogram(mapA,mapB,ax(ii,3));    % Cross-corr for the two trials.
        end
        gra_multilabel(hFig, 'row', getdatanames(data,'cell','short',cellInd));
        gra_multilabel(hFig, 'col', [getdatanames(data,'trial','long',trInd), {'Across'}]);
    case 'within'
        hFig = gra_multiplot(length(cellInd),length(cellInd),'graphborder',repmat(0.1, [1 4]));   
        ud=get(hFig, 'userdata');   ax=ud.axes_handle_array;
        for ii=1:length(cellInd)
            for jj=ii+1:length(cellInd)
                mapA=data.rate_maps(mapInd).maps{trInd,cellInd(ii)};
                mapB=data.rate_maps(mapInd).maps{trInd,cellInd(jj)};
                plotCorrelogram(mapA,mapB,ax(ii,jj));
            end
        end
        if prms.addSpatialMap
            for ii=1:length(cellInd)
                gra_plotmap(data.rate_maps( rates_mapmatch(data,SD.selMap) ).maps{trInd,cellInd(ii)},'handle',ax(ii,ii));
            end
        end
        gra_multilabel(hFig, 'row', getdatanames(data,'cell','short',cellInd));
        gra_multilabel(hFig, 'col', getdatanames(data,'cell','short',cellInd));
        gra_multilabel(hFig, 'title', {data.trials(trInd).trialname});
        
end

function plotCorrelogram(mapA,mapB,hAxis)
% Plots the correlogram in an axis.
cgram=map_crosscorr(mapA,mapB);
[gridness,grProps] = map_gridprops(cgram,'crossCorrMode',1);
cgram(isnan(cgram)) = min(min(cgram)) - ( (max(max(cgram))-min(min(cgram)))*0.1 ); % Auto-scale unvisited bins to first value.
colormap(gra_tintcolormap(10));
imagesc(cgram,'parent',hAxis);   axis(hAxis,'square');   axis(hAxis,'off');

hold(hAxis, 'on');
% Gridness Text %
fontspec = {'fontunits', 'normalized', 'fontsize', 0.12};
% Vector to closest peak %
v = [0 0; grProps.crossCorrClosestPeak];            % output form map_gridprops is relative to centre.
v = v + repmat(  fliplr(ceil(size(cgram))./2), 2, 1 );
pH=plot(hAxis,v(:,1),v(:,2),'k-');  set(pH,'linewidth',4);
% Cross-hairs %
text('string', num2str(gridness, '%3.2f'), 'position', [0 1], 'units', 'normalized','verticalalignment', 'top', 'horizontalalignment', 'left', fontspec{1:end},'color','k','parent',hAxis);
plot(hAxis, [size(cgram,2)/2 size(cgram,2)/2], get(hAxis,'ylim'), 'k--');
plot(hAxis, get(hAxis,'xlim'), [size(cgram,1)/2 size(cgram,1)/2], 'k--');
hold(hAxis, 'off');        
        
        
        
        
        
        
        