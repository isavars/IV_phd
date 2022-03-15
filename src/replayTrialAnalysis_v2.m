function [ Res, prms ] = replayTrialAnalysis_v2(varargin)
% Run the replay analysis on a trial-by-trial basis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis parameters %
prms.masterListPath = '/home/laurenz/localHD/replay_analysis/replayMasterTrialList.xlsx'; %'/home/tomw/dataLocalHDD/replay_analysis/replayMasterTrialList.xlsx';  %  'D:\tempAdultData\testData\MasterTrialList.xlsx'; % '/home/laurenz/localHD/replay_analysis/replayMasterTrialList.xlsx';
%%% (1) Parameters for this function %%%
prms.nFieldsSpkProps = 5;
% Definition of complex spikes %
prms.spkWidthThrForCS = 0.3;  % 0.3, 25 and 5 are all good values, from a histogram of the whole dataset, 2018-04-10. Best prediction is spkWidth vs ACMom, after this, rate only removes 2-3 extra cells.
prms.ACMomThrForCS    = 25;   
prms.rateThrForCS     = 5;
prms.xcZeroBinThr     = 3;
% Definition and detection of brain state: these are the frequency bands that define where theta and delta can potentially be. %
prms.thetaSearchBand = [5 10];
prms.deltaSearchBand = [2 3.5];
% Minimum number of CS cells for dataset (otherwise just don't run analysis) %
% Important that this is applied *before* rate filtering or checking for cross-tet duplicates,
% so the actual number in the final table ('nCellInDecode') will be lower than this.
prms.minCSInExpt = 27;
% Which trials to include? %
prms.analyseFirstRUNOnly = 1;
% Debugging %
prms.testDataMode    = 0;
prms.rippTestFig     = 0;
prms.brainStatePlots = 0;

%%% (2) Definition of sleep/wake states: goes to function 'getBrainStateHardThr' %%%
% Time windows for analysis %
prms.posSR       = 50;          % Sample rates: I think these are completely standard, but there is the option ..
prms.eegSR       = 250;         %    ..  to override these, if the need came about.
prms.stateWindow = 1.6;         % States calculated in sliding windows this wide ..
prms.stateStep   = 0.8;         % .. that slide by this much. 1.6 and 0.8 (sec) is the Csicsvari/Dupret default.
prms.slidingSpeedWins = 0;      % Are speed windows sliding? (can make them not, possibly fair as speed already smoothed).
prms.minEpochDur = [];         % Brain state is only defined for states that are consistent for at least this long.
% Parameters for spectrogram and frequency ratios %phaseAnalysisSquareTrack_phasePrec
prms.thetaFr     = [];          % Theta frequency to use for theta/delta ratio - can be caller-supplied, but if empty, this function calculates using FFT.
prms.thetaBW     = 1.5;         % Band width in which theta poEEGwer is calculated: centred on prms.thetaFr.
prms.deltaFr     = [];          %
prms.deltaBW     = 2;
prms.nTapers     = 2;           % 
% Thresholds for state definitions %
prms.speedThrRun   = 2.5;     % For run, must have speed faster than this ..
prms.TDRatioThrRun = 2;     %  .. and TD ratio greater than this (makes sure we are in a 'theta' epoch).
prms.speedThrREM   = 1;       % For sleep states, speed must be less than these ..
prms.speedThrSWS   = 2.5;     % .. REM needs a stricter threshold, otherwise RUN easily misclassified as REM.
prms.TDRatioThrSWS = 2;    % SWS TD ratio must be less than this
prms.TDRatioThrREM = 2;       % REM TD ratio must be greater than this.

%%% (3) Detection of ripples fromsdfsdf power envelope: goes to function 'detectRipples_react' %%%
prms.rippFiltOrder = 500;
prms.rippFiltFr1   = 100;
prms.rippFiltFr2   = 250;
prms.egfSampRate   = 4800;  % I think that this is universal, and could be completely hard-coded, but I have left it as a parameter here as a reminder that I have made this assumption.
prms.nEGFsToUse    = 1;
prms.rippPowCombMethod = 'sum';
prms.ripplePowerStdThrForCh = 4;
prms.rippThrMethod     = 'perc'; % 'std' or 'perc'
prms.peakThresh        = 99;    % Threshold for classification as ripple: either stds above mean, or %-ile.
prms.minRippleDuration = 5;    % Must cross threshold for at least this long, in ms.

%%% (4) Params for LM's 'DACQ_dataLoader' function (includes making rate maps) %%%
prms.trueRoot = '/home/laurenz/localHD/replay_analysis/rawData/'; %'/home/tomw/dataLocalHDD/replay_analysis/rawData/';   %  'D:\tempAdultData\testData\';% '/home/laurenz/localHD/replay_analysis/rawData/'  % This is the actual 'root', in the sense of the folder above the rat folders.
prms.rootDir  = '';       % This one will be fed to LM's function, it refers to the rat folder, and will therefore change each loop.
% prms for pos files
prms.loadPos = 1;
prms.ScalePosPPM = 400; %ppm
prms.posMaxSpeed = 400; %in cm/s
prms.posSmooth = 400; %in ms
prms.posHead = 0.5;
% prms for tetrode loading
prms.loadTet = 1;
prms.cutTag2 = 'fin';
prms.cutTag1 =  '_sqTrack';  %   '_sqTrack'; %   '_novLinTrack'; % 
% prms for eeg loading
prms.loadEEG = 1;
prms.downsampleEGF = 0; %y/n
% prms for map making
prms.makeMaps = 1;  % We use the open field maps to pre-screen for cross-tet recordings.
prms.mapType = 'rate'; % 'dir'; 'rate'
%rate maps
prms.placeBin    = 10; %cam pix/spatial bin
prms.smoothStyle = 'box';  %  'adapt'; % 'adapt' - adaptive smoothing; 'box' - boxcar
prms.placeSmooth = 5;  % size smoothing kernel (in bins)
prms.adSmooth    = 200;  % alpha parameter - don't change
prms.speedFilt   = [2.5 400]; % lower/uppper limits cm/s. %% This RM param is 'active' as it also applies to the linearised rate maps for decoding.
% % path scaling (will only work for standard square atm)
prms.scalePath     = 0; %y/n
% prms.minOccForEdge = 50; %in pos samples (i.e. 1s)
% prms.boxExtent     = 250; %n of cam pix for 62.5 cm box
% prms.envCheck      = {'famBox','CCEMorph','CCEsq'}; % for which envs to do path scaling
% Text progress on or off? %
prms.verbose   = 0;
%%% FINISH Params for LM's 'DAcellCQ_dataLoader' function %%%

%%% Params for MUA vector calculation - goes to 'calculate_MUA'.
prms.MUA.Zscore     = 0; %y/n
prms.MUA.kernelType = 'Gaussian'; % 'Gaussian'; 'boxcar';
prms.MUA.GaussWinSD = 10; % SD of smoothing window for MUA, in ms.
prms.MUA.boxcarWin  = 30;
prms.MUA.BinWidth   = 1; %in ms
%%% Params for MUA burst detection - goes to detectRipples_v2, using 'MUA' mode. %%%%
prms.MUA.startThresh       = 0;
prms.MUA.peakThresh        = 3;
prms.MUA.maxRippleDuration = 750;
prms.MUA.minRippleDuration = 80;
prms.MUA.minSep            = 25;
prms.MUA.cellThresh        = 0.15;
prms.MUA.absMinNCells      = 5;
prms.MUA.startBlank        = 0;

%%% Parameters for linearisation of rate maps (goes to lineariseLinearTrack or phaseAnalysisLinearisePosAndGetDirs)
prms.LinPos.trackLength  = 600;                 % This is how many pixels the track should become, *after* it has been scaled.
prms.LinPos.runDimension = 2;                   % This is the dimension to run (X=1, Y=2)
prms.LinPos.dirTolerance = circ_ang2rad( 70 );  % Tolerance for heading direction, relative to arm direction, for calculating run direction.
prms.LinPos.filtSigmaForRunDir = 1.5;  % Original sqTrack was 3s
prms.LinPos.durThrCohRun       = 5;    % Original sqTrack was 10s
prms.LinPos.durThrJump         = 1;    % Original sqTrack was 2s.

%%% linear map params
prms.LinMaps.binSize           = 10; % bin size (pixel) - will give 2.5cm bins
prms.LinMaps.smooth            = 1;  % y/n
prms.LinMaps.smoothKernelSD    = 2;  % SD, in bins
prms.LinMaps.speedFilter       = 1;  % y/n
prms.LinMaps.speedFilterLimits = prms.speedFilt; % See setting above - lock settings for linear and open field together.
prms.LinMaps.posIsCircular     = 1;
prms.LinMaps.remTrackEnds      = 0;  % Remove this many bins from each end of the linear track. Set to 0 to remove none.
prms.LinMaps.minNSpksInRateMap = 50; % min number of spikes for linear rate maps
prms.LinMaps.minPeakRate       = 0;

%%% decoding - this goes to function 'decodePos_replayAnalysis_v3' %%%
prms.decoding.replayEventDetect = 'joint';      %  'mua';  %   'ripple'; 'joint';
prms.decoding.replayType        = {'preplay','replay'};    %    {'awake','preplay','replay'};    %  
prms.decoding.mapType           = 'dir'; %'dir' - use best fit from CW or CCW runs;'all' - use full mapsprms.decoding.Tlength_SWR = 0.02; %in s 
prms.decoding.Tlength_SWR       = 0.02; %in s
prms.decoding.useOverlap        = 1;
prms.decoding.winshift          = 0.01;
prms.decoding.binSizeLin        = prms.LinMaps.binSize; %copy over for ease
prms.decoding.sBinSzMap   = 2.5; %spatial bin size in cm
prms.decoding.binRange    = 6;            % yRange and ..
prms.decoding.spdBnd      = [0, 2500];    %  .. spdBnd for CB lineTraj function.
prms.decoding.doShuffle   = 1;
prms.decoding.nShuff      = 100; %
prms.decoding.shuffleType = 'postTime';     %    % 'cellID'; %  'postTime';  %    'map';  %    'postPos';  %
prms.decoding.cutTag1     = prms.cutTag1;  % Need to copy 'cutTag1' into the 'decoding' sub-struct, as it tells DC function whether pos is circular (_sqTrack) or linear (_novLinTrack).
prms.decoding.Tonline     = 0.3;

% Some special settings:
prms.GDReplication        = 1;   % Refers to all map, event, decoding etc settings, not the shuffling, which is dealt with separately.
% These two go together:
prms.crossMatchPseudoData = 1;
prms.decoding.runOnlineDC = 0;

% Some specific settings that always need to change if analysing linear track: these
% are re-set here, so to go from 'standard' SqTr to LT only requires one param change.
if strcmp( prms.cutTag1, '_novLinTrack' )
    prms.LinMaps.posIsCircular = 0;
    prms.LinPos.durThrCohRun   = 2;   % Coherent run dur shorter (can run end-to-end track faster)
end

% GD replication settings.
if prms.GDReplication
    prms.decoding.replayEventDetect = 'mua';
    prms.decoding.Tlength_SWR       = 0.02; %in s
    prms.decoding.useOverlap        = 0;
    prms.decoding.binRange          = 2;
    prms.MUA.peakThresh             = 2;
    prms.MUA.maxRippleDuration      = 800;
    prms.MUA.minRippleDuration      = 100;
    prms.MUA.GaussWinSD             = 15; % SD of smoothing window for MUA, in ms.
    prms.LinMaps.binSize            = 4; % 10=2.5cm bins. Matches bin N, not bin size. 4=1cm bins
    prms.LinMaps.smoothKernelSD     = 5;  % SD, in bins
    prms.LinMaps.speedFilterLimits  = [5 400];
    prms.LinMaps.minNSpksInRateMap  = 10; % min number of spikes for linear rate maps
    prms.LinMaps.minPeakRate        = 1;
    prms.LinPos.durThrCohRun        = 2.5; % Duraction of coherent runs in one direction - not specified by GD, CB paper 5s, taking a guess here.
    prms.speedThrSWS                = 2;   % Used when defining brain state
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);
            dotInd = strfind( varargin{ii}, '.' );
            if isempty(dotInd)
                prms.(varargin{ii}) = varargin{ii+1};
            else
                prms.( varargin{ii}(1:(dotInd-1)) ).( varargin{ii}((dotInd+1):end) ) = varargin{ii+1};
            end
        end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1a. Load in the trial list, create a list of individual experiemnts.
if prms.testDataMode
    prms.masterListPath = 'replayMasterTrialList.xlsx';   % Load the local version from in the code folder.
end
TrialList        = readtable( prms.masterListPath );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% % QUICK HACK, 27/07/18 - run only new datasets (collected June 2018 onwards).
if 0
    TrialList = TrialList(428:end,:);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


TrialList.ratNum = cellfun( @(x) str2double(x(2:end)), TrialList.rat_ID );               % A numerical rat identifier is also helpful - need to remove the 'r'.
TrialList        = TrialList(    strcmp( TrialList.cutFile_ID1, prms.cutTag1 ),    : );    % Remove the reactivation expts, just leave the sqTrack expts.

if prms.testDataMode
    dataToKeepInd = false( size( TrialList, 1 ), 1 );
    dataToKeepInd = dataToKeepInd | (TrialList.ratNum==579 & TrialList.age==21);
    TrialList     = TrialList( dataToKeepInd, : );
    prms.trueRoot = 'E:\Dropbox\matlab\testDataReact\';   %  C:\Users\Tom Wills
end
[~,~,expUIDInd] = unique( [TrialList.ratNum, TrialList.experiment_ID], 'rows' );  % 'expUIDInd' is a unique identifier for experiemnts (as the field 'experiemnt_ID' in the sheet is only unique within animal).
uniqueExpList   = unique(expUIDInd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1b. Pre-alocate results table
if prms.analyseFirstRUNOnly
    maxNTrial = 3;
else
    maxNTrial = max(  accumarray( expUIDInd, ones(size(expUIDInd)) )  );
end
cellDum   = cell(  1, maxNTrial     );
if 1     % For neatness, table field name declarations are hidden in this folded IF
    varList   =   {
             'rat',         nan; ...   
             'age',         nan; ...
             'nCellInDecode',   nan; ...
 
             'trial',      cellDum; ...    % One value per trial, but this is a cell string array.
             'trialDur',   nan( 1, maxNTrial ); ...
                          
             'spkWidth',    cell(1,1); ...  % These fields will all hold a (1,nCell) vector, 1 per expt 
             'acMom'        cell(1,1); ...  % (spike props are mean values).
             'MRbyExpt',    cell(1,1); ...
             'cellNum',     cell(1,1); ...
             'tet' ,        cell(1,1); ...
             'xcZeroBin',   cell(1,1); ...
             'cellPairID',  cell(1,1); ...
             'mapCorr',     cell(1,1); ...
             'spdHistRun',  cell(1,1); ...
             
             'MRInRunEpochByCell', cell( 1, 1 ); ...
             'MRInRunEpoch',       nan( 1, 1 );  ...
             'MRInRunEpochSEM',    nan( 1, 1 );  ...
             'MRInSWSEpochByCell', cell(1, maxNTrial); ...
             'MRInSWSEpoch',       nan( 1, maxNTrial ); ...
             'MRInSWSEpochSEM',    nan( 1, maxNTrial ); ...
             
             
             'isSleepTrial',    nan( 1, maxNTrial ); ...
             'nEGFs',           nan( 1, maxNTrial ); ...
             'meanEGFSD',       nan( 1, maxNTrial ); ...
             'sumEnv99Prc',     nan( 1, maxNTrial ); ...
             'swrPowerAbs',     nan( 1, maxNTrial ); ...
             'thetaFreq',       nan( 1, maxNTrial ); ...
             'deltaFreq',       nan( 1, maxNTrial ); ...
             'dirSim',          nan( 1, maxNTrial ); ...
             'speed',           cell(1, maxNTrial); ...
             'meanSpdRun',         nan( 1, 1 ); ...
             'meanSpdRunNoStops',  nan( 1, 1 ); ...
             'medSpdRun',          nan( 1, 1 ); ...
             'medSpdRunNoStops',   nan( 1, 1 ); ...
%              'meanMUA',         nan( 1, maxNTrial ); ...
%              'thrMUA',          nan( 1, maxNTrial ); ...
             
             'eventTimesSWR',   cellDum; ...
             'eventTimesMUA',   cellDum; ...
             'eventTimesJoint',   cellDum; ...
             'freqSpectSWR',   cellDum; ...
             'stateTimes',     cellDum; ...
%              'rawSWR',   cellDum; ...
%              'filtSWR',   cellDum; ...           
             
             'spikeTimes',   cellDum; ...
              
             'lin_rMaps',   cellDum; ...
             'lin_pMaps',   cellDum; ...
             'linPos',      cellDum; ...
             'linDir',      cellDum; ...
             'brainStates',   cellDum; ...

            
             % The following fields are created in decodePos_replayAnalysis, but I think they need to also be declared here %
             'AwDCMedErr', 	    nan( 1, maxNTrial ); ...                                                                    %
             'AwDC75thErr',     nan( 1, maxNTrial ); ...
             
             'actPosAtEvent', cell( 1, 1 ); ...
             
             'bstRes',      cellDum; ....
             'bstSpd',      cellDum; ...
             'bstYInterc',  cellDum; ...
             'wCorr',       cellDum; ...
             'jDist',       cellDum; ...
             'jStd',        cellDum; ...
             'jMean',       cellDum; ...
             
             
             'bstResShuf',  cellDum; ...
             'bstSpdShuf',  cellDum; ...
             'bstYIntercShuf',  cellDum; ...
             'wCorrShuf',   cellDum; ...
             'jDistShuf',   cellDum; ...
             'jStdShuf',    cellDum; ...
             'jMeanShuf',    cellDum; ...
             
             'pVals',              cellDum; ...
             'decodeSpikeTimes',   cellDum; ...

                 };
end
dummyTableRow                          = cell2table( varList(:,2)' );   % We can't just create the whole table here - we need a dummy row
dummyTableRow.Properties.VariableNames = varList(:,1);
Res                                    = repmat( dummyTableRow, length(uniqueExpList), 1 );
Res.Properties.UserData                = prms;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Load the data, and run pre-processing necessary for replay decode: 
%    find SWS brain state, ripples, MUA events, make rate maps.
tic
parfor ii=1:length(uniqueExpList)   % ii=iterator for EXPERIMENT 
    
        TrialListForExp = TrialList(  expUIDInd==uniqueExpList(ii),  :   );

        %%% QUICK HACK INSERTION: if you want to just run one dataset, uncomment this line 
        %%%                       and insert the correct rat number and age (adult age = 35).
%         if ~any(TrialListForExp.age(1)==[30]) || ~any(TrialListForExp.ratNum(1)==[739]);  continue;   end
        
%         try
            ResForExp       = loadAndPreProcessSubFunc( TrialListForExp, dummyTableRow, prms );
            Res(  ii,  :  ) = ResForExp;
%         catch ME
%             fprintf(1, 'Error with dataset r%d P%d \n', TrialListForExp.ratNum(1), TrialListForExp.age(1) );
%             errME{ii,1} = ME;
%         end

        
end
Res = Res( ~isnan(Res.rat), : );
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2b) [Post-bublication insertion]. If requested, carry out cross-rat RUN-sleep
%    pseudo-data control. The function below re-matches sleep spike/event time data 
%    from one session with rate maps from another, age matched, always across rats, 
%    maps or spikes randomly sub-sampled to match N.
if prms.crossMatchPseudoData
    Res = Res( Res.nCellInDecode>=25, : ); % Get canonical dataset first (don't want to use low cell sessions in cross-match).
    Res = replayMakePseudoData( Res, prms );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Run the event-by-event decode. Each dataset run serially, the parfor is applied
%    at the event level for loop inside decodePos_replayAnalysis_v3.
rng default;  % IMPORTANT: make sure we always get the same set of shuffles. 
tic
trToDecode = false(1,3);
for ii=1:length(prms.decoding.replayType);   trToDecode = trToDecode | strcmp(prms.decoding.replayType{ii},{'preplay','awake','replay'});   end
% errDS = {};
for ii=1:size(Res,1)
    
    nEvent = sum( cellfun( @length, Res.eventTimesJoint(ii,trToDecode) ) );
    fprintf( 1, 'Decoding dataset %d out of %d (%d events). \n', ii, size(Res,1), nEvent );
    
%     try
        ResDecode = decodePos_replayAnalysis_v3( Res( ii, : ), prms.decoding );
        v = ResDecode.Properties.VariableNames;
        for jj=1:length(v)
            Res.(v{jj})( ii, : )  = ResDecode.(v{jj});
        end
%     catch
%         errDS{end+1} = [Res.rat(ii) Res.age(ii)];
%     end
     
end
toc

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub-function that loads the data, finds brain state, finds event times %
function [ResForExp] = loadAndPreProcessSubFunc( TrialList, dummyTableRow, prms )
% This sub-function loads the data for one whole experiemnt, and then loops through trials to generate 
% within-trial cell-pair correlations.
% Note that the input 'TrialList' (table of trial information derived from excel master sheet) contains
% only the rows (row=trial) relevant to this experiment.

% Are we going to analyse all RUN trials, or only the first?
% (Note - if we ever did want to analyse all, would probably be better
% to treat each RUN as a separate dataset, as the sleep in the middle
% would otherwise have two sets of results, replay against RUN1 and 
% preplay against RUN2, which would be confusing to deal with).
if prms.analyseFirstRUNOnly
    [~, trialTypes] = cellfun( @(x) strtok(x,'_'), TrialList.trialName, 'UniformOutput', 0 );
    RUNTrialInd     = ~strcmp( trialTypes, '_sleepHP' );    %    strncmp( trialTypes, '_sqTrack', 8 ) | strcmp( trialTypes, '_novlinTrack' );   TW: changed this May 2019 to load open field FAM-NOV expts.
    firstRunInd     = find( RUNTrialInd , 1, 'first' );     % NOTE: adult env is called _sqTrackLarge, so compare first 8 characters only.
    TrialList       = TrialList(  [-1 0 1]+firstRunInd,   :   );
    nTrial          = 3;
else
    nTrial    = size(TrialList,1);
end
ResForExp = dummyTableRow;   % Need to build table rows independently for expts, due to PARFOR.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1a. First up, we need to load in all of the cell spiking properties, for all trials in the experiment,
%     so as we can decide  which are complex spike, and thereby know which tets have complex spikes recorded on them.
path           = [prms.trueRoot, TrialList.rat_ID{ 1 }, filesep];
spkPropsStruct = load( [path, strtok(TrialList.trialName{1}, '_'), TrialList.cutFile_ID1{1}, '_spikeProps.mat'] );
for jj=2:nTrial
    spkPropsStruct(jj) = load( [path, strtok(TrialList.trialName{ jj }, '_'), TrialList.cutFile_ID1{ jj }, '_spikeProps.mat'] ); % The naming convention is (date)(trial sequence letter)(cut file ID1[as in crib sheet])(_spikeProps), e.g. 170302d_sqTrack_spikeProps
end
spkPropsArray = cat(3, spkPropsStruct(1:end).spkProps );
spkPropsMean  = nanmean( spkPropsArray, 3 );
isCSCellInd   = spkPropsMean(:,3)>=prms.spkWidthThrForCS & spkPropsMean(:,4)<=prms.ACMomThrForCS & spkPropsMean(:,5)<=prms.rateThrForCS;
isCSTetList   = unique( spkPropsMean( isCSCellInd, 2 ) );
% IMPORTANT: if there aren't enough CS cells, quit already at this point %
if sum( isCSCellInd ) < prms.minCSInExpt;
    fprintf( 1, 'Rat %s P%d: only %d CS cells, quitting \n', TrialList.rat_ID{ 1 }, TrialList.age( 1 ), sum( isCSCellInd ) );
    return;   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1c. Load in the raw data - for all trial in this experiment %
prms.rootDir = [prms.trueRoot, TrialList.rat_ID{ 1 }, filesep];
prms.Tnames  = TrialList.trialName( : );
%     prms.Tets    = 1:8; %%% COMMENTED OUT BY LM - WILL NOW AUTOMATICALLY FIND ALL EXiSTING TETRODES
[ Pos, Spk, EEG, Maps ,dataLog] = DACQ_dataLoader( 'hc', prms ); % LM function.
% index for CS in Spk struct - TW note to LM: why did you put this 'isCS' variable in? I assumed you could directly read across cell identity from the spkProps to the tint cuts,
% but it seems like 'isCS' has been put there in the assumption that you can't. Is that correct?
[~,~,isCS] = intersect( spkPropsArray(isCSCellInd,1:2,1), Spk(1).CellTetN,'rows' );
isCS       = sort(isCS);
% filter out cells that don't fire enough in run trial (i.e. won't have rate map)
runInd     = cellfun('isempty',strfind( TrialList.trialName, 'sleepHP' ));
spikeTimeN = cellfun(@length,Spk(runInd).spikeTimes(isCS)) > prms.LinMaps.minNSpksInRateMap;
isCS       = isCS(spikeTimeN);
% Final filtering out of non CS from all trials %
for jj=1:length(Spk) 
    Spk(jj).spikeTimes = Spk(jj).spikeTimes(isCS);
    Spk(jj).CellTetN   = Spk(jj).CellTetN(isCS,:);
end
Maps.rMaps   = Maps.rMaps( :, isCS );
spkPropsMean = spkPropsMean( isCS, : );
% Calculate 0-bin x-correlogram of cells, not on same tetrode, looking for cross-recorded cells %
cnt = 1;
[ResForExp.xcZeroBin{1}, ResForExp.mapCorr{1}] = deal(  nan(   (length(Spk(runInd).spikeTimes)^2 - length(Spk(runInd).spikeTimes))/2,  1 )   );
ResForExp.cellPairID{1} = nan(   (length(Spk(runInd).spikeTimes)^2 - length(Spk(runInd).spikeTimes))/2,  2 );
for kk=1:length(Spk(runInd).spikeTimes)
    for mm=(kk+1):length(Spk(runInd).spikeTimes)
        if map_spatialcorr( Maps.rMaps{runInd, kk}, Maps.rMaps{runInd, mm} ) < 0.8;   continue;   end
        xcTemp                         = spk_crosscorr( Spk(runInd).spikeTimes{kk}, Spk(runInd).spikeTimes{mm}, 0.001, 0.05, length(Pos(runInd).XY)/prms.posSR );
        ResForExp.xcZeroBin{1}(cnt,1)  = xcTemp(51) ./ mean( xcTemp( [1:49 52:101] ) );
        ResForExp.cellPairID{1}(cnt,:) = [kk mm];
        ResForExp.mapCorr{1}(cnt,1)    = map_spatialcorr( Maps.rMaps{runInd, kk}, Maps.rMaps{runInd, mm} );
        cnt                            = cnt+1;
    end
end
% Remove cells defined as cross-tetrode recorded - this is slightly involved, as we need to 
% be careful we don't remove more cells than necessary, so go down the list and remove those 
% cells that appear in the most pairs.
xRecPair = ResForExp.cellPairID{1}( ResForExp.xcZeroBin{1} > prms.xcZeroBinThr, : );
cellToRemList = [];
for ii=1:size(xRecPair,1)
   if isnan(xRecPair(ii,1));   continue;   end
   if sum( xRecPair(:)==xRecPair(ii,1) )  >=  sum( xRecPair(:)==xRecPair(ii,2) )
       cellToRem = xRecPair(ii,1);
   else
       cellToRem = xRecPair(ii,2); 
   end   
   xRecPair( any(xRecPair==cellToRem,2), : ) = nan;
   cellToRemList = [cellToRemList, cellToRem];
end
xRecFiltInd = true( size( Spk(1).spikeTimes ) );
xRecFiltInd( cellToRemList ) = false;
for jj=1:length(Spk) 
    Spk(jj).spikeTimes = Spk(jj).spikeTimes(xRecFiltInd);
    Spk(jj).CellTetN   = Spk(jj).CellTetN(xRecFiltInd,:);
end
spkPropsMean           = spkPropsMean( xRecFiltInd, : );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1c. For each run trial, get single, best, low-sample rate theta trace: needed for theta cycle definition and theta/delta ratio %
%     For sleep trials, get the delta frequency.
%     This section needs to go before the main trial loop, because we really want to estimate delta from sleep and theta from run,
%     *before* calculating theta/delta ratio for each single trial in the main analysis.
[thetaFr, deltaFr] = deal( nan(1,nTrial) );
bestEEG            = cell( 1, nTrial );
for jj=1:nTrial
    % Is trial jj sleep or run? Need to set frequency and speed filter parameters accordingly %
    if strcmp( TrialList.trialName{jj}(end-6:end), 'sleepHP' )
        freqBand    = prms.deltaSearchBand;
        speedFilter = find( Pos(jj).speed < 2.5 );
    else
        freqBand    = prms.thetaSearchBand;
        speedFilter = find( Pos(jj).speed > 2.5 );
    end
    speedFilterEEG     = reshape(   bsxfun(@plus, speedFilter.*5, -4:0)',  [],  1  );  % Need to be careful that what is reshaped has the right row-col orientation: the filter index needs to come out with following indices contiguous!
    % For each trial, loop through all the EEG channels, to find the best SNR in the required range. %
    [peakFr, SNR] = deal( nan(length(EEG(jj).eegData),1) );
    for kk=1:length(EEG(jj).eegData)
        if ~isempty( EEG(jj).eegData{kk} ) && any( EEG(jj).tetEEG(kk)==isCSTetList )
            [peakFr(kk), ~, SNR(kk)] = eeg_powerspec( EEG(jj).eegData{kk}(speedFilterEEG), prms.eegSR, 'hfCutOff', 25, 'thetaBand', freqBand, 'fig', 0, 'freqRes', 0.02);
        end
    end
    % Find the channel with best SNR, record the frequency of interest and save the actual best EEG trace for later use %
    [~,bestEEGInd] = nanmax(SNR);
    bestEEG{jj}    = EEG(jj).eegData{ bestEEGInd };
    if ~strcmp( TrialList.trialName{jj}(end-6:end), 'sleepHP' )
        thetaFr(jj)  = peakFr( bestEEGInd );
    else
        deltaFr(jj)  = peakFr( bestEEGInd );
    end
end
% Record these estimates for later cross-checking %
ResForExp.thetaFreq(1,1:length(thetaFr)) = thetaFr;
ResForExp.deltaFreq(1,1:length(deltaFr)) = deltaFr;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Now make the linearised rate maps. Rate maps are based on linearised position,
%    and always come in two sets: filtered for CW and CCW runs.
RUNTrInd = 2; % (Note, if we ever wanted to analyse more than one RUN trial per dataset, we would need to loop from this point onward.)
if strcmp( prms.cutTag1, '_sqTrack' )
    linPos  = phaseAnalysisLinearisePosAndGetDirs( Pos(RUNTrInd).XYraw, Pos(RUNTrInd).direction, prms.LinPos );
elseif strcmp( prms.cutTag1, '_novLinTrack' )
    linPos  = lineariseLinearTrack( Pos(RUNTrInd).XYraw, Pos(RUNTrInd).direction, prms.LinPos );
end
dirInd                             = [linPos.cohRunInd == 1, linPos.cohRunInd == 2];
[lin_rMaps,lin_pMaps,RMCellsInUse] = makeLin_rMaps_v3(linPos.linPos, dirInd, Spk(RUNTrInd).spikeTimes, Pos(RUNTrInd).sRate, Pos(RUNTrInd).speed, prms.LinMaps);
ResForExp.lin_rMaps{1,RUNTrInd}    = lin_rMaps;
ResForExp.lin_pMaps{1,RUNTrInd}    = lin_pMaps;
% Record 'dirSim', overall correlation between CW and CCW rate maps %
ResForExp.dirSim(1,RUNTrInd)       = map_spatialcorr( lin_rMaps{2}, lin_rMaps{3} );
% These are needed for online decoding, at the next step %       
ResForExp.linPos{1,RUNTrInd}       = linPos.linPos;
ResForExp.linDir{1,RUNTrInd}       = linPos.cohRunInd;
% 2b. Filter spike times again, so we retain only those cells contributing to rate maps.
for jj=1:length( Spk )
    Spk(jj).spikeTimes = Spk(jj).spikeTimes( RMCellsInUse );
end
spkPropsMean       = spkPropsMean( RMCellsInUse, : );
% 2c. Filter the *spikes* by coherent run epoch, to get a measure of mean firing rate when running.
% First get a list of run epoch times.
for ii=1:2
    runEpochs{ii,1} = find( diff([0; linPos.cohRunInd==ii; 0]) == 1 ) - 1;  % Run starts for dir ii
    runEpochs{ii,2} = find( diff([0; linPos.cohRunInd==ii; 0]) == -1 );     % Run ends for dir ii
end
runEpochs = cell2mat( runEpochs );                  % Get a unified list of CW and CCW epochs ..
runEpochs = sortrows( runEpochs, 1 );               % .. sorted in time.
runEpochs = runEpochs .* (1/Pos(RUNTrInd).sRate);   % Convert to seconds.
% .. now filter the spike times for each cell by run epoch times, and store mean rate in epochs.
MRInRunEpoch = zeros(length( Spk(RUNTrInd).spikeTimes ), 1);
for itCl=1:length( Spk(RUNTrInd).spikeTimes )
    spkGTEpoch          = bsxfun( @gt, Spk(RUNTrInd).spikeTimes{itCl}, runEpochs(:,1)' );
    spkLTEpoch          = bsxfun( @lt, Spk(RUNTrInd).spikeTimes{itCl}, runEpochs(:,2)' );
    spkInEpochInd       = any(  spkGTEpoch & spkLTEpoch,   2  );
    MRInRunEpoch(itCl) = sum( spkInEpochInd ) ./ sum( diff(runEpochs,1,2) );   
end
ResForExp.MRInRunEpochByCell{1}  = MRInRunEpoch;
ResForExp.MRInRunEpoch(1)        = mean( MRInRunEpoch );
ResForExp.MRInRunEpochSEM(1)     = std( MRInRunEpoch ) / sqrt( length(MRInRunEpoch) );




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2c. Get metadata *for the experiment*. This needs to come after all filtering/spike exclusion,
%     so the final numbers and cell IDs match with what is in the replay decode.
ResForExp.rat(1,1)           = TrialList.ratNum( 1 );
ResForExp.age(1,1)           = TrialList.age( 2 );
ResForExp.trial(1,1:nTrial)  = TrialList.trialName;
ResForExp.nCellInDecode(1,1) = length( Spk(1).spikeTimes );
ResForExp.spkWidth{1,1}      = spkPropsMean(:,3);
ResForExp.acMom{1,1}         = spkPropsMean(:,4);
ResForExp.MRbyExpt{1,1}      = spkPropsMean(:,5);
ResForExp.cellNum{1,1}       = spkPropsMean(:,1);
ResForExp.tet{1,1}           = spkPropsMean(:,2);
ResForExp.meanSpdRun(1)         = nanmean( Pos(runInd).speed ); 
ResForExp.meanSpdRunNoStops(1)  = nanmean( Pos(runInd).speed(  Pos(runInd).speed>=2.5  ) );
ResForExp.medSpdRun(1)          = nanmedian( Pos(runInd).speed ); 
ResForExp.medSpdRunNoStops(1)   = nanmedian( Pos(runInd).speed(  Pos(runInd).speed>=2.5  ) );
ResForExp.spdHistRun{1}         = histc( Pos(runInd).speed, 0:2.5:50 );
% Assign spike times to output. %
for jj=1:length( Spk )
    ResForExp.spikeTimes{1,jj} = Spk(jj).spikeTimes;
    ResForExp.trialDur(1,jj)   = dataLog(jj).trialDur;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Now loop over trials (PRE, RUN, POST) and detect SWR and MUA burst periods.
%    These will define the 'events' that will be decoded.
for jj=1:nTrial  % jj=iterator for TRIAL. TW: changed this May 2019 to load open field FAM-NOV expts. Was '1:3'.
    
    ResForExp.speed{1,jj} = Pos(jj).speed; % store raw speed 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3a. Define SWS times, based on movement and theta/delta ratio %
    %     Assuming sleep/run trials alternate, we get the frequency type for this trial type from this trial,
    %     and the frequency of the other band from the mean of that taken from the preceding and following trials.
    prms.thetaFr = nanmean(   thetaFr( max([1 jj-1]) : min([nTrial jj+1]) )   );
    prms.deltaFr = nanmean(   deltaFr( max([1 jj-1]) : min([nTrial jj+1]) )   );
    allStates    = getBrainStateHardThr( Pos(jj).speed, bestEEG{jj}, prms );
    ResForExp.brainStates{1,jj} = zeros( size(allStates.sws) );
    ResForExp.brainStates{1,jj}( allStates.sws ) = 1;
    ResForExp.brainStates{1,jj}( allStates.rem ) = 2;
    ResForExp.brainStates{1,jj}( allStates.run ) = 3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3b. Find ripples in SWS.
    % Filter the EGF traces in the ripple band, and convert to a trace of power in this band. %
    % First, to make sure we only use the EGFs we want, discard those that are (a) filtered LP+Notch,
    % b) higher gain (there are often 'low gain' and 'high gain' duplicate sets of EEGs - we want only the low).
    % NOTE - this has all now gone in a sub-function.
    [eventTimesSWR, ~, rippPeakPowers, EGFProps] = detectRipples_react( EEG(jj), isCSTetList, prms );  % NOTE - eventTimes in sec at this point.
    % Record some characterisitics of EGFs %
    ResForExp.nEGFs(1,jj)       = EGFProps.nEGFs;
    ResForExp.meanEGFSD(1,jj)   = EGFProps.meanEGFSD;
    ResForExp.sumEnv99Prc(1,jj) = EGFProps.sumEnv99Prc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3c. Make sure that we are only using windows from the correct brain state (SWS for sleep)
    % Get the correct state vector for the trial. This is a logical index with timebase prms.stateStep: also need to convert to a (nEpoch,2) array of state epoch start and end times.
    stateInd = allStates.sws;
    stateDiff  = diff( [0 double(stateInd) 0] );
    stateStart = find( stateDiff==1 );
    stateEnd   = find( stateDiff==-1 ) - 1;
    stateTimes = [stateStart', stateEnd'] .* prms.stateStep;
    stateTimes = bsxfun( @plus, stateTimes, [-1 1].*(prms.stateStep/2) );  % Push the state starts and ends 'out' by half a state bin, so they run between step edges not step centres.
    % Now, test which spiking windows (windowTimes) lie within a valid state epoch (epochTimes), and discard those which aren't.
    winStartInState = bsxfun( @gt, eventTimesSWR(:,1), stateTimes(:,1)' );
    winEndInState   = bsxfun( @lt, eventTimesSWR(:,2), stateTimes(:,2)' );
    winInStateEpoch = winStartInState & winEndInState;          % This is a (nWin, nStateEpoch) logical array, true where a win lies within a state epoch.
    validWinInd     = any(  winInStateEpoch ,   2   );          % However, we are only interested in whether win is valid (not which epoch), so just do ANY along rows.
    eventTimesSWR   = eventTimesSWR( validWinInd, : );
    % Assign output for SWRs, also get some analysis diagnostics for ripple detection %
    ResForExp.swrPowerAbs(1,jj)   = mean( rippPeakPowers( validWinInd ) );
    ResForExp.eventTimesSWR{1,jj} = eventTimesSWR;
    ResForExp.freqSpectSWR{1,jj}  = EGFProps.rippSpect( :, validWinInd );
    ResForExp.stateTimes{1,jj}    = stateTimes;
%     ResForExp.rawSWR{1,jj}        = EGFProps.allRippRaw( :, validWinInd );
%     ResForExp.filtSWR{1,jj}       = EGFProps.allRippFilt( :, validWinInd );
    % Filter *spikes* by SWS epochs, so that we can get a measure of the Mean Rate specifically in SWS.
    MRInSWSEpoch = zeros(length( Spk(RUNTrInd).spikeTimes ), 1);
    for itCl=1:length( Spk(RUNTrInd).spikeTimes )
        spkGTEpoch         = bsxfun( @gt, Spk(RUNTrInd).spikeTimes{itCl}, stateTimes(:,1)' );
        spkLTEpoch         = bsxfun( @lt, Spk(RUNTrInd).spikeTimes{itCl}, stateTimes(:,2)' );
        spkInEpochInd      = any(  spkGTEpoch & spkLTEpoch,   2  );
        MRInSWSEpoch(itCl) = sum( spkInEpochInd ) ./ sum( diff(stateTimes,1,2) );   
    end
    ResForExp.MRInSWSEpochByCell{1,jj} = MRInSWSEpoch;
    ResForExp.MRInSWSEpoch(1,jj)       = mean( MRInSWSEpoch );
    ResForExp.MRInSWSEpochSEM(1,jj)    = std( MRInSWSEpoch ) / sqrt( length(MRInSWSEpoch) );
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3d. Define MUA bursts %
    % First need to create a filter of valid times (resolution=1ms) to look for MUAs, based on brain state %
    MUAStateFilter = false( 1, dataLog(jj).trialDur*1000 );  
    stateTimesAsInd = ceil( stateTimes .* 1000 );
    stateTimesAsInd(stateTimesAsInd==0) = 1;
    for kk=1:size(stateTimes,1)
        MUAStateFilter(   stateTimesAsInd(kk,1) : stateTimesAsInd(kk,2)   ) = true;
    end
    % Now calculate MUA and define burst-events %
    [ MUA, ~, MUA_nSpikes ]            = calculate_MUA( Spk(jj).spikeTimes, dataLog(jj).trialDur, 1000, prms.MUA);
    prms.MUA.nCells                    = length(Spk(jj).spikeTimes);  % These two are needed by the ..
    prms.MUA.MUA_nSpikes               = MUA_nSpikes;                 % detectRipples_v2 function, in MUA mode.
    [eventTimesMUA,~,~,meanMUA,thrMUA] = detectRipples_v2( MUA, 'MUA', 1000, MUAStateFilter, [], prms.MUA );  % From this function, eventTimesMUA is in units of ms (integers - can be used as index))
    ResForExp.eventTimesMUA{1,jj}      = eventTimesMUA;
%     ResForExp.meanMUA(1,jj)            = meanMUA;
%     ResForExp.thrMUA(1,jj)             = thrMUA;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3e. Also define the MUA+SWR coincidence events - specifically, the MUA event windows which overlap in any part with an SWR window  %
    eventTimesSWR                   = ceil( eventTimesSWR .* 1000 );  % Convert real secs to integer ms, for compatibility with 'joinRipples' function.
    ResForExp.eventTimesJoint{1,jj} = joinRipples( eventTimesSWR, eventTimesMUA, 'BA'); %join ripples - 'BA' mode means we get the *times* based on MUA, whenever one of these *overlaps* at all with an SWR.

end  % END for Sleep Trials

end

    












