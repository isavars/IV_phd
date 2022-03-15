function []=gra_mapfig_series()
% Draw maps based on a timeseries of data from within trial.

SD=gss;
data=evalin('base',SD.selData{1});

% Establish if we are doing multiple cells for one trial, or multiple trials with one cell %
if length(SD.selTrial)==1
    dataInd=1:length(SD.selCell);
    rowLabels=getdatanames(data,'cell','short',SD.selCell);
    figTitle=[SD.selData(1) {data.trials(SD.selTrial(1)).trialname}];
elseif length(SD.selCell)==1
    dataInd=1:length(SD.selTrial);
    rowLabels=getdatanames(data,'trial','short',SD.selTrial);
    figTitle=[SD.selData(1) getdatanames(data,'cell','long',SD.selCell)];
else
    error('Please select either one trial or one cell - you CAN select multiple cells in one trial, or multiple trials for one cell, but NOT multiple cells in multiple trials');
end

% Get the user defined parameters for the map series %
prms=rates_params;

% If user asked for Pos maps, then make sure that this will be shown, regardless of if there are multiple cells selected %
if strcmp(prms(1).mode,'pos') && length(SD.selTrial)==1
    dataInd = 1;
    rowLabels=getdatanames(data,'trial','short',SD.selTrial);
    rowLabels{1} = {rowLabels{1}, 'Positions'};
end
    
hFig=gra_multiplot(length(dataInd),length(prms)+1,'figborder',[0.5 0.25 0.5 0.5]);
ud=get(hFig,'userdata'); hAx=ud.axes_handle_array;

for ii=1:length(prms)+1
    %%% Make the maps %%%
    if ii==1
        prmsTemp=prms(1);  prmsTemp.filt_time=[];
        maps=rates_main(data,prmsTemp);
        colLabels{1}='All';
    else
        maps=rates_main(data,prms(ii-1));
        colLabels{ii}=[num2str(prms(ii-1).filt_time(1)) '-' num2str(prms(ii-1).filt_time(2))];
    end   
    %%% Plot %%%
    if strcmp(prms(1).mode,'pos')
        cellIndex = 1;
    else
        cellIndex = SD.selCell;
    end
    for jj=dataInd
        if length(SD.selTrial)==1
            gra_plotmap(maps{ SD.selTrial(1), cellIndex(jj) },'handle',hAx( jj, ii ));
        elseif length(SD.selCell)==1
            gra_plotmap(maps{ SD.selTrial(jj), cellIndex(1) },'handle',hAx( jj, ii ));
        end           
    end
end
gra_multilabel(hFig, 'row', rowLabels);
gra_multilabel(hFig, 'col', colLabels);
gra_multilabel(hFig, 'title', figTitle);