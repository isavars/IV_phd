function []=gra_distributivehypothesis(SD)
% Plot directional rate map and matching distributive hypothesis map.
% Shows log difference ratio in corner.

ratioShift=0.01;  % Add this number to both maps when calculating DH mean ratio.

% Make maps %
data=evalin('base',SD.selData{1});
params=rates_params('bin',50,'smooth',1,'bin_dir',15,'smooth_dir',1);
params.space='dir';
params.alg='pd';
maps=rates_main(data,params);
params.alg='dh';
mapsDh=rates_main(data,params);

% Figure %
hFig = gra_multiplot(length(SD.selCell), length(SD.selTrial), 'graphborder', repmat(0, [1 4]), 'figborder', [1 0.5 1 0.5],'plotsize',[2 2]);
ud=get(hFig,'userdata');
set(hFig,'defaulttextfontname','arial','inverthardcopy','off','color','white','numbertitle','off','name',['Dist Hypo: ' SD.selData{1}]);
uimenu('label', 'Save as eps ..', 'callback', 'scan_base(''saveEps'', gcf);');

% Plot %
for ii=1:length(SD.selTrial)
    for jj=1:length(SD.selCell)
            mapDh=mapsDh{SD.selTrial(ii),SD.selCell(jj)};
            map=maps{SD.selTrial(ii),SD.selCell(jj)};
            mapDh=mapDh.*( mean(map)./mean(mapDh) );      % Mean normalise maps
            
            axes(ud.axes_handle_array(jj,ii));
            hMap=gra_plotmap(mapDh,'handle',ud.axes_handle_array(jj,ii),'text_pos','none');
            set(hMap,'color','y');
            hold on
            gra_plotmap(map,'handle',ud.axes_handle_array(jj,ii));
            hold off
            
            % Ratio %
            dhr=mean( abs( log( (maps{ii,jj}+ratioShift) ./ (mapDh+ratioShift) ) ));
            text('position', [0 0.85], 'units', 'normalized', 'HorizontalAlignment', 'left', 'string', num2str(dhr,'%4.2f'), ...
                'FontUnits', 'normalized', 'VerticalAlignment','cap', 'fontsize', 0.10, 'color',[0 0 0]);
            
    end
end
    
    
    
    
    
    
    
    
    
    
