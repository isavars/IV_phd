function []=gra_clusterviewer(SD)
% View clusters and tetrodes of currectly selected cells.

hFig=gra_multiplot(2,3,'figborder',ones(1,4).*0.1,'graphborder',zeros(1,4));
ud=get(hFig, 'userdata');   hAxArray = ud.axes_handle_array;

chKey = {[1 2],[1 3],[1 4]; [2 3],[2 4],[3 4]}; % Key for which tetrode channels to plot in each axes - R,C layout corresponds to axes layout.
clusterCMap = [0.7 0.7 0.7; 0 0 1; 0 1 0; 1 0 0; 1 0 1; 0 1 1; 0 0.5 0; 0.7 0.7 0; 0.4 0.4 0.4];    % Tint cluster colours. Only have first 8 so far.

data=evalin('base',SD.selData{1});
selTet=data.trials( SD.selTrial(1) ).cells( SD.selCell(1) ).tet;
cellIndex=cat(1,data.trials(1).cells.tet) == selTet;
tetDataTemp=data.trials( SD.selTrial(1) ).cells( cellIndex );

% Sort cells by number. This has the effect of moving cell 0 to the front, so it plots first, and therefore at the back. %
[~,cellSortInd] = sort( cat(1,tetDataTemp.cellnum) );
tetDataTemp = tetDataTemp(cellSortInd);

for ii=1:2
    for jj=1:3
        for kk=1:length(tetDataTemp)
            if isempty(tetDataTemp(kk).st);   continue;   end
            xPlot = tetDataTemp(kk).wf_amps(:, chKey{ii,jj}(2) );
            yPlot = tetDataTemp(kk).wf_amps(:, chKey{ii,jj}(1) );
            if tetDataTemp(kk).cellnum == 0; m=2; else m=2; end
            line('xdata',xPlot,'ydata',yPlot,'linestyle','none','marker','s','markersize',m,'markerfacecolor',clusterCMap( tetDataTemp(kk).cellnum+1, : ),'markeredgecolor','none','parent',hAxArray(ii,jj));
        end
        set(hAxArray(ii,jj),'xlim',[0 255],'ylim',[0 255],'box','on','xtick',[],'ytick',[]);
    end
end

set(hFig,'name',[SD.selData{1} ' - ' data.trials(SD.selTrial(1)).trialname ' - Tet ' num2str(selTet)]);