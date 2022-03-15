
function [H] = gra_wf(varargin)

% Draw all mean maximum waveforms in a data set.
%
%       gra_wf('plotType');         (Call from SCAn GUI);
%       gra_wf(data, 'plotType');   (Call from command line);
%
% plotType - 'allchannels', 'maxchannel'
 
SD=gss;

if nargin==1
    SD = gss;   data = evalin('base', SD.selData{1});
    plotType = varargin{1};
else
    data = varargin{1};
    plotType = varargin{2};
end

n_trial = length(SD.selTrial);
n_cell = length(SD.selCell);

if strcmp(plotType, 'maxchannel');
    n_figs = 1; 
    n_row = n_cell;
    n_col = n_trial;
elseif strcmp(plotType, 'allchannels');
    n_figs = 1:length(SD.selCell);
    n_row = 4;
    n_col = n_trial;
end


%%% Make the figure(s) %%%
for kk=n_figs

	H = gra_multiplot(n_row, n_col, 'axesborder', [0.1 0.1 0.1 0.1], 'figborder', [1 0.5 1 0.5]);
	set(H,'paperorientation', 'landscape', 'defaulttextfontname', 'arial');
	set(H,'inverthardcopy','off','color','white');      % Workaround for Matlab bug that prints white text as black when 'axis off' is used.
    fig_ud = get(H, 'userdata');
    
    
	%%% Draw plots into axes %%%
    switch plotType
        case 'maxchannel'
            for ii=1:n_col   
                for jj=1:n_row
                    cTemp=data.trials(SD.selTrial(ii)).cells(SD.selCell(jj));
                    gra_plotwf(cTemp.wf_means,cTemp.scalemax,'max',fig_ud.axes_handle_array(jj,ii));
                end
            end        
            % Labels for trial and cell %   (This needs to go after drawMap, or text disappears. Don't know why).
            gra_multilabel(H, 'row', getdatanames(data, 'cell', 'short', SD.selCell));
            gra_multilabel(H, 'col', getdatanames(data,'trial','long', SD.selTrial));

        case 'allchannels'
            %%%% All channels %%%%
            for ii = 1:n_col     
                ah = fig_ud.axes_handle_array(:,ii);
                cTemp=data.trials(SD.selTrial(ii)).cells(SD.selCell(kk));
                gra_plotwf(cTemp.wf_means, cTemp.scalemax, 'all', ah);
            end 
            %%% Labels for trial and cell %%%   (This needs to go after drawMap, or text disappears. Don't know why).
            gra_multilabel(H, 'title', ['Cell ', num2str(data.trials(1).cells(SD.selCell(kk)).cellnum), ' (t', num2str(data.trials(1).cells(SD.selCell(kk)).tet), ')']);
            gra_multilabel(H, 'row', {'Ch1', 'Ch2', 'Ch3', 'Ch4'});
            gra_multilabel(H, 'col', getdatanames(data,'trial','long',SD.selTrial));
    end
        

end





