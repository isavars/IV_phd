
function [counts] = spk_isihist(spkTr1, spkTr2, maxInt, varargin)

% Calsulates the inter-spike interval histogram, not the cross-correlogram.
%
%       [hist] = spk_isihist(spktr1, spktr2, maxInterval);
%
% Max interval is in ms.
%
% Or with options ...
%
%       [rtn] = spk_isihist(spktr1, spktr2, maxInterval, 'paramName', paramValue);
%
% 'trialDur', trialLength   - Trial length, in Sec. Y Axis will be spikes/s,
%                             rather than counts.
% 'divisions', N            - N divisions in histogram. Default = 50. (For
%                             one side of histogram only).
% 'hAxis', axisHandle       - Pass an axis handle to draw into.
% 'yLabel', 1 or 0          - Draw or not Y axis labels (default: draw)
% 'plot', 1 or  0           - Draw a figure;


% Input Options %
opt.hAxis = [];
opt.trialDur = [];
opt.divisions = 50;
opt.yLabel = 1;
opt.drawPlot = 1;
for ii=1:2:length(varargin)
    opt.(varargin{ii}) = varargin{ii+1};
end

% Work in ms %
spkTr1 = spkTr1.*1000; 
spkTr2 = spkTr2.*1000;

% Find all Intervals (N-(N+1), N-(N+2), N-(N+3) etc .. ) %
nSpk = length(spkTr1);
int = repmat(0,[nSpk nSpk-1]);
for ii=1:nSpk-1                         % Missing out last iteration prevents 0 intervals (N-N).
    ind = [nSpk-ii+1:nSpk 1:nSpk-ii];   % Index to shift spikes to N+nSpk position (wrapping around first-last spikes)
    spkTr2Shift = spkTr2(ind);
    intTemp = spkTr1 - spkTr2Shift;
    int(:,ii) = intTemp;
end       

% Bin %
binwidth = maxInt/opt.divisions;
bins = -maxInt:binwidth:maxInt;
if isempty(spkTr1) || isempty(spkTr2)
    % This is necessary for graceful failure when nSpk=0 %
    counts = repmat(NaN, size(bins));
else
    % This is normal %
    counts = histc(reshape(int,1,[]),bins);
end
if ~isempty(opt.trialDur)
    counts = (counts ./ opt.trialDur);
end

% Plot %
if opt.drawPlot
    if ~isempty(opt.hAxis)
        axes(opt.hAxis)
    else
        figure;
    end
    hBar = bar(bins,counts,'histc');
    set(hBar, 'edgecolor', 'none', 'facecolor', [0 0 0]);
    set(gca, 'xlim', [0 maxInt], 'xtick', 0:maxInt/10:maxInt, 'xticklabel', ...
        {'0','','','','',num2str(maxInt/2),'','','','',num2str(maxInt)}, ...
        'fontunits', 'normalized', 'fontsize', 0.1, 'tickdir', 'out');
    if ~opt.yLabel
        set(gca, 'ytick', []);
    end
end

% {num2str(-maxInt),'','','','',num2str(-maxInt/2),'','','','','0','','',''
% ,'',num2str(maxInt/2),'','','','',num2str(maxInt)}



