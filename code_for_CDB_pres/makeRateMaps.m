function [rMaps, sm_spkMaps, sm_pMaps, sACs] = makeRateMaps(spkTimes, positions, sampleTimes, ppm, varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% params
prms.smooth       = 'adaptive'; % 'boxcar; 'adaptive'
prms.smoothKernel = ones(5,5);
prms.binSize      = 2.5; % in cm^2
prms.Fs           = 50;
prms.alpha        = 200;

prms.makeSACs = 0;

if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end   
% in case just 1 cell
if ~iscell(spkTimes)
    spkTimes = {spkTimes};
end

%% posMap
positions = ceil( bsxfun(@minus,positions,min(positions,[],1)) + eps ); % re-order to pos range to start from 1 in both dimensions
binSizePix = round(ppm/100 * prms.binSize); % this many pix in one bin
posBinned = ceil(positions ./ binSizePix);
posMapRaw = accumarray(posBinned(~isnan(posBinned(:,1)),:),1,[nanmax(posBinned(:,1)) nanmax(posBinned(:,2))]) ./ prms.Fs;
unVisPos = posMapRaw == 0; % keep record of invisited positions

%% make maps

rMaps = cell(length(spkTimes),1);
sm_spkMaps = cell(length(spkTimes),1);
sACs = cell(length(spkTimes),1);
if strcmp(prms.smooth,'boxcar')
    sm_pMaps = {};
else
    sm_pMaps = cell(length(spkTimes),1);
end

for i = 1:length(spkTimes)
    % spike Map
    spkPosBinInd = arrayfun(@(x) find(sampleTimes - x > 0,1,'first'),spkTimes{i},'UniformOutput',0); % as sample times can be somewhat irregular we can't just bin by sample rate
    spkPosBinned = posBinned(cell2mat(spkPosBinInd),:);
    spkMapRaw    = accumarray(spkPosBinned(~isnan(spkPosBinned(:,1)),:),1,[max(posBinned(:,1)) max(posBinned(:,2))]);
    
    % smooth
    switch prms.smooth
        case 'boxcar'
            % smooth
            if ~isempty(sm_pMaps)
                sm_pMaps{1} = filter2(prms.smoothKernel,posMapRaw);
            end
            sm_spkMaps{i} = filter2(prms.smoothKernel,spkMapRaw);
            %rate map
            rMaps{i} = sm_spkMaps{i} ./ sm_pMaps{1};
            rMaps{i}(unVisPos) = NaN;
            
        case 'adaptive'
            [rMaps{i}, sm_spkMaps{i}, sm_pMaps{i}] = rates_adaptivesmooth(posMapRaw,spkMapRaw,prms.alpha); % SCAN function
            
        otherwise
            error([prms.smooth ' is not a valid option for smoothing rate maos']);
    end
    
    if prms.makeSACs
        sACs{i} = map_crosscorr(rMaps{i}, rMaps{i});    
    end
    
end

% 
% subplot(1,3,1);
% imagesc(rMap,[0 nanmax(rMap(:))]); colormap(jet); axis square
% subplot(1,3,2);
% spk_crosscorr(spkTimes,'AC',0.005,0.5,1300,'plot',gca); axis square
% 
% [spatCorr] = map_crosscorr(rMap, rMap);
% subplot(1,3,3);
% imagesc(spatCorr); colormap(jet); axis square
end

