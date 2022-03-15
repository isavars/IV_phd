 function behavFilterPlot( varargin )
% This function will read in muiltiple DACQ data files and then plot those
% periods where the animal is showing some specific behaviour (e.g. 
% immobility).
%
% Note: immobility is only filter implemented at the moment
%
%  Usage:  behavFilterDACQ
%          behavFilterDACQ( optionalInputStruct )
%          behavFilterDACQ( 'inputName', inputVal, .. etc .. )
%
%  Optional inputs/analysis parameters (supply as " ,'fieldname',value, " comma-separated list) :
%
%          prms.behavFilter = 'immobility'; % 'immobility';
%          %for pos loading
%          prms.ScalePosPPM = 400; % scale to this pix/m;
%          prms.posHead = 0.5; % pos of head with respect to LEDs
%          prms.posMaxSpeed = 400; % upper speed limity in cm/s
%          prms.posSmooth = 0.4; % in s
%          prms.smoothSpeedSD = 1; % kernel SD for smoothing in s - needs to be integer!
%          prms.maxSpeed = 1.5; % speed max for immobility (cm/s)
%          prms.briefMovementT = 0.5; % in s; periods of this length that comtain movement will be ignores
%          prms.minTimeL = 10; % in s; min duration for behavioural periods  
% 
%          prms.debug = 1; %y/n - will plot a figure with overview over periods that were accepted /rejected
%
%  LM 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% params
prms.fNameOut = 'r661_180226mno_sleepHP';
prms.behavFilter = 'immobility'; % 'immobility';
%for pos loading
prms.ScalePosPPM = 400; % scale to this pix/m; 
prms.posHead = 0.5; % pos of head with respect to LEDs
prms.posMaxSpeed = 400; % upper speed limity in cm/s
prms.posSmooth = 0.4; % in s
prms.smoothSpeedSD = 1; % kernel SD for smoothing in s - needs to be integer!
prms.maxSpeed = 1.5; % speed max for immobility (cm/s)
prms.briefMovementT = 0.5; % in s; periods of this length that comtain movement will be ignores
prms.minTimeL = 10; % in s; min duration for behavioural periods 

prms.savePlot = 0;

%parse optional inputs
if nargin>0 && isstruct(varargin{1})
    optIn = varargin{1};
    f=fieldnames(optIn);
    for i=1:length(f);   prms.(f{i}) = optIn.(f{i});   end
elseif nargin>0 && ischar(varargin{1})
    for i=1:2:length(varargin)
        prms.(varargin{i}) = varargin{i+1};
    end
end

% open new figure 
figure;

%% pre stuff
% directory
defaultDir = 'Z:\IsabellaV\recording_data\'; %'D:\tempAdultData\'; %'Z:\lmuessig\!postDoc\recording_data\CA1\'; 'D:\tempAdultData\';
if ~isdir(defaultDir)
    defaultDir = 'C:\';
end

% grab set files
[setNames, dataDir] = uigetfile([defaultDir '*.set'],'Select .set Files of Trials for Which to Filter Data by Behaviour','MultiSelect','on'); 
if ischar(setNames)
    setNames = cellstr(setNames); % in case just one trial is selected
end
% fail gracefully
if ~ischar(dataDir)
    warning('Loading aborted! Neieieieieieiein...');
    return
end
%% get relavant periods in trials
for i = 1:length(setNames)
    %load set - taken from SCAN
    [sFile, sFileText] = load_set([dataDir setNames{i}(1:end-4)]); %load set file
    % in case there is overhang in data
    trialDur = sFile.duration;
    if ~~rem(trialDur,10)
        trialDur = trialDur - rem(trialDur,10);
    end
    switch prms.behavFilter
        case 'immobility'
            %load pos 
            pFile = load_pos([dataDir setNames{i}(1:end-4)],sFile.tracked_spots);
            pFile.led_pos = floor(pFile.led_pos .* (prms.ScalePosPPM/pFile.pixels_per_metre)); % scale pix/m
            pFile.pixels_per_metre = prms.ScalePosPPM;
            pFile.header{strmatch('pixels_per_metre',pFile.header(:,1)),2} = num2str(prms.ScalePosPPM);  % For POSTPROCESS_POS_DATA
            %Interpolate, Filter and Smooth %
            [~, ~, speed] = postprocess_pos_data(pFile, prms.posMaxSpeed, prms.posSmooth, sFileText, prms.posHead);
            % smooth speed with Gaussian
            kernelBins = prms.smoothSpeedSD * pFile.sample_rate;
            kernel = fspecial('Gaussian',[1 6*kernelBins], kernelBins); % kernel
            speedSM = imfilter(speed', kernel, 'conv','replicate'); % smooth
            % this will yield speed averages of length 'prms.smoothSpeedSD' 
            % windows overlapping by length 'prms.smoothSpeedSD/2'
            stepSz = prms.smoothSpeedSD/2 * pFile.sample_rate;
            speedSM = repmat(speedSM(1:stepSz:end),stepSz,1); % upsample to pos sample rate
            speedSM = speedSM(:); 
            % find periods of immobility
            logicSpeedFilt = speedSM <= prms.maxSpeed;
            startInds = find(diff(logicSpeedFilt) > 0) + 1; % start points
            endInds = find(diff(logicSpeedFilt) < 0); % end points
            % check if first or last pair is incomplete
            if endInds(1) < startInds(1)
                startInds = [1; startInds];
            end
            if startInds(end) > endInds(end)
                endInds = [endInds; length(speedSM)];
            end
            % start/end times in seconds - rounded inclusively
            timeStart = floor(startInds / pFile.sample_rate); 
            timeStop =  ceil(endInds / pFile.sample_rate); 
            if startInds(1) == 1
                timeStart(1,1) = 0;
            end
            % in case there is overhang, remove it
            if timeStop(end) > trialDur
                timeStop(end) = trialDur;
                % for plotting also remove from speed vector
                speedSM = speedSM(1:trialDur*pFile.sample_rate);
                speed = speed(1:trialDur*pFile.sample_rate);
            end
            % exclude periods of brief movements
            briefMovInd = timeStart(2:end) - timeStop(1:end-1) < prms.briefMovementT;
            timeStart(find(briefMovInd)+1) = []; timeStop(find(briefMovInd)) = [];
            % remove periods which are too short
            durInd = timeStop - timeStart < prms.minTimeL;
            timeStart(durInd) = [];
            timeStop(durInd) = [];
            
            % plot periods that were extracted from each trial
            logicSpeedFiltPlot = false(length(speedSM),1);
            if ~isempty(timeStart)
                startPlot = timeStart * pFile.sample_rate; endPlot = timeStop * pFile.sample_rate;
                
                if startPlot(1) == 0; startPlot(1) = 1; end
                
                for f = 1:size(startPlot,1)
                    logicSpeedFiltPlot(startPlot(f):endPlot(f)) = true;
                end
            end
            xVals = 1:length(speedSM);
            subplot(length(setNames),1,i);
            hold on
            plot(speed,'k-','linewidth',0.5);
            tempX = xVals; tempY = speedSM;
            tempX(logicSpeedFiltPlot) = NaN; tempY(logicSpeedFiltPlot) = NaN;
            plot(tempX,tempY,'r-','linewidth',2);
            tempX = xVals; tempY = speedSM;
            tempX(~logicSpeedFiltPlot) = NaN; tempY(~logicSpeedFiltPlot) = NaN;
            plot(tempX,tempY,'g-','linewidth',2);
            xAxTix = 0:pFile.sample_rate*60:length(speedSM);
            xAxStr = cellstr(num2str((xAxTix/pFile.sample_rate)'));
            set(gca,'XTick',xAxTix,'XTicklabel',xAxStr,'xlim',[0 length(speedSM)],'ylim',[0 40],'ytick',[0 20 40],'yticklabel',{'0','20','40'});
            hold off
            title([setNames{i} '; Duration of period: ' num2str(sum(diff([timeStart timeStop],[],2)))],'interpreter','none');
            if i == length(setNames)
                xlabel('Time (s)');
            end
            clear endPlot startPlot logicSpeedFiltPlot tempX tempY xAxStr xAxTix        
    end
end
if prms.savePlot
    print(gcf,'-painters','-bestfit','-loose','-dpdf',[dataDir prms.fNameOut '.pdf']);
end
end


