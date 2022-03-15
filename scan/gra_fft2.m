function []=gra_fft2(input1,varargin)
% Plot rate maps and 2-D fourier transforms
%
% To call whilst running scanning GUI:
%
%       gra_fft2(gss,shufFFTMax);   or
%       gra_fft2(gss,[]);           if shuffled data FFT max populations are not available (no thresholding of FFT will be done).
%
% shufFFTMax should have the structure shufFFTMax.shuf{1:nDataset,2}, where column 1 is a cell array of N Datasets worth of 3-d arrays,
% each with the format (1:nTrial,1:nCell,1:nShuffles), and column 2 is a cell string array, 1:nDataset, with the dataset name.
%
% Calling this function from code is not working at the moment ... (20-03-2012)


% Parse input %
if isfield(input1,'selData')
    % If called by SCAn GUI %
    datasetName=input1.selData{1};
    data = evalin('base',datasetName);
    trInd=input1.selTrial;
    cellInd=input1.selCell;
    rateMapsForDisplay = data.rate_maps( rates_mapmatch(data,input1.selMap) ).maps;
    binSize=data.rate_maps( rates_mapmatch(data,input1.selMap) ).params.bin;
    if isempty(rates_mapmatch(data,'bin',binSize,'smooth',1))
        rateMapsForFFT = rates_main(data,rates_params('bin',binSize,'smooth',1));
    else
        rateMapsForFFT = data.rate_maps( rates_mapmatch(data,'bin',binSize,'smooth',1) ).maps;
    end
    shufFFTMax=varargin{1};
else
    % If called in code %
%     data = input1;
%     trInd=varargin{1};
%     cellInd=varargin{2};
%     rateMaps = varargin{3};
%     binSize=varargin{4};
%     shufFFTMax=varargin{5};
end

hFig=gra_multiplot(length(cellInd),length(trInd)*3);     hAx=getappdata(hFig,'axesHandles');

fftPad = 256;
powerSpecSigThr = 95;   % Refers to the percentile of max power specs from spatially shuffeld data.

if ~isempty(shufFFTMax)
    % If there is shuffled data, take the following parameters from this %
    angBinSize= shufFFTMax.params.angBinSize;                
    normMode = shufFFTMax.params.normMode;             
    lowFreqLimit = shufFFTMax.params.lowFreqLimit;       
    highFreqLimit = shufFFTMax.params.highFreqLimit;
else
    % If there is none, use these as defaults.
    angBinSize= 3 ;                % Angular bin for radial mean power
    normMode = 'peak';             % Rate normalisation
    lowFreqLimit=2.5;       % units are cycles per metre.   1.67 cyc/m => 60cm WL,    2.5cyc/m => 40cm WL (1.5 cyc per box length)
    highFreqLimit=6.45;     %  ..                           8    cyc/m => 12.5cm WL, 6.45cyc/m => 15.5cm WL (4 cyc per box length)
end
    
% Check that FFT and rate map parameters match for shuffled and real data %
if ~isempty(shufFFTMax)
    if shufFFTMax.params.binSize~=binSize
        disp('WARNING! Shuffled data was not created using the same bin size as those in the real data you have selected - you shouldn''t trust the signifcance threshold.');
    end
end

trial_count = 1;
for ii=trInd 
    cell_count = 1;
    % Set frequecy limits for this trial %
    binSpatFreq = 1 / ((1/data.trials(ii).ppm) * binSize);
    fftLimitOuter =  fftPad / (binSpatFreq/highFreqLimit);
    fftLimitInner =  fftPad / (binSpatFreq/lowFreqLimit);
    [x, y]=meshgrid(-(fftPad/2)+1:fftPad/2, -(fftPad/2)+1:fftPad/2);
    fftFreqMask = (x.^2+y.^2)>=fftLimitInner^2 & (x.^2+y.^2)<=fftLimitOuter^2;
    for jj=cellInd
        
        % Plot Rate Map %
        gra_plotmap(rateMapsForDisplay{ii,jj},'handle',hAx(cell_count,trial_count));
        
        
        %%%% Make and plot 2-D FFT %%%%%
        map=rateMapsForFFT{ii,jj};
        % If normMode~='', make some attempy to control for overall firing rate of cell %
        switch normMode
            case 'peak'
                map = map./nanmax(map(:));
            case 'mean'
                map = map./nanmean(map(:));
            case 'parseval'
            	map = map./ sqrt( nansum( map(:).^2 ));
        end
        map=map - (mean(map(~isnan(map))));     % Mean normalise
        map(isnan(map))=0;                      % Set non-visited to zero 
        ps=abs(fft2(map, fftPad, fftPad));      % zero padding the rate map
        ps=fftshift(ps);                        % shifting low frequency components to the middle
        psMax = nanmax(ps(fftFreqMask));            % Get maximum power in relevant frequency band
        % Set non-significant spectral values to 0 %
        sigThr=0;
        psThr=ps;
        if ~isempty(shufFFTMax) % shufFFTMax can be passed as [], as a placeholder, if calling without shuffled FFT data.
            dsInd=strcmp(shufFFTMax.shuf(:,2),datasetName);
            sigThr = prctile(squeeze(shufFFTMax.shuf{dsInd,1}(ii,jj,:)),powerSpecSigThr);
            if psMax >= sigThr
                psThr(ps < sigThr) = 0;    % Remove non-sig components (only in sig maps, otherwise leave all).
            end
        end
        if psMax >= sigThr;  lineSpec='k:';   else   lineSpec='r-';   end   % Change line style for plots below depending on if sig or not
        % Crop power spectrum for display (to Upper frequency limit) %
        cen=fftPad/2;
        cr=round(fftLimitOuter);
        psCrop = psThr( cen-cr+1:cen+cr, cen-cr+1:cen+cr );        % 
        % Plot FFT %
        fftAx=hAx(cell_count,trial_count+1);
        imagesc(psCrop,'parent',fftAx);     % Plot
        hold(fftAx,'on');
        % Plot Cross hairs %
        crCen=size(psCrop,1)/2;     % Centre of cropped ps
        plot(fftAx,[1 size(psCrop,2)],(ones(1,2).*crCen)+1,lineSpec);  % Cross-hairs. Note that the centre of the FFT comes
        plot(fftAx,(ones(1,2).*crCen)+1,[1 size(psCrop,1)],lineSpec);  % out 1 pixel right and down of the centre of the matrix.
        % Plot Inner and outer frequency-of-interest limit circles %
        th=linspace(0,2*pi,360);
        r=ones(size(th)).*fftLimitOuter;  [x,y]=pol2cart(th,r);  x=x+crCen+1;  y=y+crCen+1;  plot(fftAx,x,y,lineSpec);
        r=ones(size(th)).*fftLimitInner;  [x,y]=pol2cart(th,r);  x=x+crCen+1;  y=y+crCen+1;  plot(fftAx,x,y,lineSpec);
        % Print values for max FFT and for FFT significance threshold %
        fontspec = {'fontunits','normalized','fontsize',0.1,'units','normalized','backgroundcolor','w','parent',fftAx};
        text('string',num2str(round(psMax),'%d'),'position',[0 1],'verticalalignment','top','horizontalalignment','left',fontspec{1:end});
        text('string',num2str(round(sigThr),'%d'),'position',[1 0],'verticalalignment','bottom','horizontalalignment','right',fontspec{1:end});
        axis(fftAx,'off')
        
        
        
        %%%%% Make and plot Angular mean of powers %%%%%%
        ps(~fftFreqMask) = nan;
        %  .. Radialise FFT ..
        [y,x]=ind2sub(size(ps),1:numel(ps));        % Get x,y coords for indices of ps.
        x=x-fftPad/2;   y=y-fftPad/2;                               % Set x=0,y=0 at centre of array
        [psAng,~]=cart2pol(x,y);                                    % Convert to radial array
        psRad=ps(1:numel(ps));                      % Linearise FFT with same index.
        % .. take angular mean and get max. %
        psAngBinned = ceil(psAng./     ((angBinSize/180)*pi)    );
        psAngBinned = psAngBinned - min(psAngBinned) + 1;
        psRadBinned = accumarray( psAngBinned', psRad, [], @nanmean );
        % Plot .. %
        angAx=hAx(cell_count,trial_count+2);
        % .. change the line appearance, depending of significance %
        angThr = prctile(squeeze(shufFFTMax.shufAng{dsInd,1}(ii,jj,:)),powerSpecSigThr);
        if nanmax(psRadBinned) >= angThr;    lineSpec='b-';   else   lineSpec='r:';   end
        angBinList = angBinSize:angBinSize:360;
        plot(angAx,angBinList,psRadBinned,lineSpec);
        % Print values for max FFT and for FFT significance threshold %
        hold(angAx,'on');
        plot(angAx,[angBinList(1) angBinList(end)], ones(1,2).*angThr, 'r-');
        set(angAx,'xlim',[0 180],'xtick',[0 180],'ytick',get(angAx,'ylim'));
        
        
        % Bump counter %
        cell_count=cell_count+1;
    end
    trial_count=trial_count+3;
end









