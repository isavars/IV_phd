function []=gra_spatialinfo(SD)
% Plot rate map and print matching spatial information in corner


% Make maps %
data=evalin('base',SD.selData{1});
params=rates_params;
if ~isempty(params.adaptive_smooth)
    [maps posMaps]=rates_main(data,params);
else
    maps=rates_main(data,params);
    params.mode='pos';
    posMaps=rates_main(data,params);
end

% Figure %
hFig = gra_multiplot(length(SD.selCell), length(SD.selTrial), 'graphborder', repmat(0, [1 4]), 'figborder', [1 0.5 1 0.5],'plotsize',[2 2]);
ud=get(hFig,'userdata');
set(hFig,'defaulttextfontname','arial','inverthardcopy','off','color','white','numbertitle','off','name',['Spatial Info: ' SD.selData{1}]);
uimenu('label', 'Save as eps ..', 'callback', 'scan_base(''saveEps'', gcf);');

% Plot %
for ii=1:length(SD.selTrial)
    for jj=1:length(SD.selCell)
            
            % Plot Map % 
            axes(ud.axes_handle_array(jj,ii));
            gra_plotmap(  maps{SD.selTrial(ii),SD.selCell(jj)}  ,'handle',ud.axes_handle_array(jj,ii));
            
            % SI %
            si=map_skaggsinfo(  maps{SD.selTrial(ii),SD.selCell(jj)}  ,   posMaps{SD.selTrial(ii)}  );
            text('position', [1 1], 'units', 'normalized', 'HorizontalAlignment', 'right', 'string', num2str(si,'%4.2f'), ...
                'FontUnits', 'normalized', 'VerticalAlignment','cap', 'fontsize', 0.10, 'color',[0 0 0]);
            
    end
end
    
    
    
    
    
    
    
    
    
    
