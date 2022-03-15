function []=gra_dirprops(varargin)
%
%

% Check for calling mode, set relevant variables appropriately %
if ~isempty(varargin)
    % Called on command line %
    data = varargin{1};
    %map_ind = varargin{2};
    opt.trials = varargin{3};
    opt.cells = varargin{4};
    fig_namer = 'Directional Props';
    mapParams =  rates_params('space','dir','bin',6,'smooth',5);
    posMapParams = rates_params('space','dir','mode','pos','bin',6,'smooth',5);
else
    % Called by SCAN UI %
    SD = gss;
    data = evalin('base', SD.selData{1});
    opt.trials = SD.selTrial;
    opt.cells = SD.selCell;
    mapParams = data.rate_maps(SD.selMap).params;
    posMapParams = mapParams;
    posMapParams.mode = 'pos';
    fig_namer = [SD.selData{1}, ' Dir Props'];    
end
% [ILMaps1 ILMaps2] = rates_interleavedmaps(data,60,'space','dir','bin',6,'smooth',5);
dirMaps = rates_main(data, mapParams);
posMaps = rates_main(data, posMapParams);


% Set up figure %
h_fig = gra_multiplot(length(opt.cells),length(opt.trials),'graphborder',[0 3.5 0 0],'figborder',[1 0.5 1 0.5]);
set(h_fig,'defaulttextfontname','arial','numbertitle','off','name',fig_namer,'inverthardcopy','off','color','white');
ud = get(h_fig, 'userdata');
% Menus %
uimenu('label', 'Save as eps ..', 'callback', 'scan_base(''saveEps'', gcf);');

%%% Draw maps %%%
for ii = 1:length(opt.trials)        
    for jj = 1:length(opt.cells)
        map = dirMaps{opt.trials(ii),opt.cells(jj)};
        posSampSpk = ceil(data.trials(opt.trials(ii)).cells(opt.cells(jj)).st.*data.trials(opt.trials(ii)).sample_rate);
        %% Draw maps into selected axis %%
        axes(ud.axes_handle_array(jj, ii));
        gra_plotmap(map);
        %%% Directional Stats %%%
        % Half-height arc %
        dirStat{1} = round( (sum( map>(max(map)/2) ) / (length(map))) * 360 );
        % Angular SD %
        dirStat{2} = nan;%circ_std( ang2rad( double(data.trials(opt.trials(ii)).dir(posSampSpk))) );
        % Spatial Info %
        [dirStat{3} dirStat{4}] = map_skaggsinfo(map,posMaps{opt.trials(ii),1});
        % KL Divergence %
        dirs = (1:length(map))';   % Any uniform vector is fine, I think.
        probMap = map./sum(map);
        probUni = repmat(1/length(map), size(map));
        dirStat{5} = dir_kldivergence(dirs,probMap,probUni);
        % Watson U Sq %
        dirStat{6} = dir_watsonusq( data.trials(opt.trials(ii)).dir, data.trials(opt.trials(ii)).dir(posSampSpk) );
        % Rayleigh Vector %
        dirStat{7} = dir_rayleighvector(map);
        % Interleaved Half-trial map stability %
%         dirStat{7} = map_spatialcorr(ILMaps1{opt.trials(ii),opt.cells(jj)},ILMaps2{opt.trials(ii),opt.cells(jj)});
        %%% Plot directional Stats on Axes %%%
        IDStr = {'HH=','SD=','SI=','SIR=','KL=','WU=','RV='};
        for kk=1:length(dirStat)
            if strcmp(IDStr{kk},'HH=')
                labelStr{kk}=[IDStr{kk} num2str(dirStat{kk})];
            else
                labelStr{kk}=[IDStr{kk} num2str(dirStat{kk},'%6.2f')];
            end
        end
        text('units','normalized','position',[1,0.5],'string',labelStr,'horizontalalignment','left','fontunits','normalized','fontsize',0.125);        
    end
end
gra_multilabel(h_fig, 'row', getdatanames(data, 'cell', 'short', opt.cells));
gra_multilabel(h_fig, 'col', getdatanames(data, 'trial', 'short', opt.trials));