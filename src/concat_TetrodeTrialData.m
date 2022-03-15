function  [dirOut, fNameOut, firingsOut] = concat_TetrodeTrialData( tetNum, mode, varargin)
% This function will concatonate tetrode data across multiple trials into 1
% joint file. We assume that trial data names are of format:
% dateStamp(american format) - trial index (letter) - underscore - .....
% The output file will then be called 'dateStamp - prms.addStr - .(tetNum)
% tint mode: will create a tetrode file like dacQ would generate
% kilosort mode: 
% mountain sort mode: 
%
%       [...] = concat_TetrodeTrialData(tetNum, mode)
%       [...] = concat_TetrodeTrialData(tetNum, mode, optionalInputStruct)
%       [...] = concat_TetrodeTrialData(tetNum, mode, 'inputName', inputVal, .. etc .. )
%
%
% Inputs:   tetNum  - tetrode number
%           mode    - 'tint' - concatonate data in TInt format;
%                     'kilosort' - generate data for clustering with
%                     kilosort (needs a bit more checking etc.)
%                     'mountainsort' - generate data for clustering with
%                     mountainsort (needs a bit more checking etc.)
%
% Optional inputs/analysis parameters (supply as " ,'fieldname',value, " comma-separated list) :
%
%           prms.fNames = names of files to be concatonated (all except
%                         extension)
%           prms.dataRootDir = directory where data for tetrodes is stored
%           prms.writeDir = directory where data will be written
%           prms.addStr = ''; string to be added to file (dateStamp(XX).(tetNum))
%           prms.leaveOverhang = 1;
%           prms.sRate = 48000; sample rate of spikes - usually 48kHz
%           prms.noiseMax = 5; (kilosort only); max amplitude for 'fake'
%                           baseline variation
%
% Outputs:  
%           fOut - file name of concatonated tetrode data  
%           dirOut - full path to dir of concatonated tetrode data 
%
% LM 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% params
prms.fNames = {'211210a_famBox','211210b_famBox','211210c_CCE','211210d_famBox','211210ef_sleepHP'}; %trials to load
prms.dataRootDir = 'Z:\IsabellaV\recording_data\Adults\r933'; %'D:\tempRecordingData\r566\'; 
prms.writeDir = 'Z:\IsabellaV\recording_data\Adults\r933';  % 'D:\'; 'D:\SandBox\dataShare\'
prms.addStr = ''; % string to be added to file (dateStamp(XX).(tetNum))
% just for tint
prms.leaveOverhang = 1; % flag for removing/leaving trial overhang
% 
prms.sRate = 48000; %sample rate
prms.addNoise = 1;
prms.noiseMax = 5; %'noise' ampl (min + max) for baseline for kilosort data

prms.wfSign = -1;
%parse optional inputs
if nargin>2 && isstruct(varargin{1})
    optIn = varargin{1};
    f=fieldnames(optIn);
    for i=1:length(f);   prms.(f{i}) = optIn.(f{i});   end
elseif nargin>2 && ischar(varargin{1})
    for i=1:2:length(varargin)
        prms.(varargin{i}) = varargin{i+1};
    end
end
% assuming beginning of filename is YYMMDD
dateStamp = prms.fNames{1}(1:6);

%% read tet data
wavFin = nan(0,50,4);
TSFin = nan(0,4);
duration = 0;
for nTrials = 1:length(prms.fNames)
    %CB code
    [header, timestamp, waveforms] = read_tetrode_file([prms.dataRootDir prms.fNames{nTrials} '.' num2str(tetNum)]);
    
    % grab trial duration    
    durInd = strcmp(header(:,1),'duration');
    tempDur = sscanf((header{durInd,2}),'%d');
    % if the overhang in DACQ files is meant to be cut need to remove
    % those samples
    if prms.leaveOverhang || ~rem(tempDur,10)
        wavFin = cat(1,wavFin,waveforms);
        TSFin = cat(1,TSFin,timestamp+duration);
    else
        tempDur = tempDur - rem(tempDur,10);
        ind = timestamp(:,1) <= tempDur;
        wavFin = cat(1,wavFin,waveforms(ind,:,:));
        TSFin = cat(1,TSFin,timestamp(ind,:) + duration);
    end
    %make sure time stamps increase continuously
    duration = duration + tempDur;    
end
clear durInd timestamp waveforms tempDur


%% write data to disk
switch mode
    
    case 'tint'
        fNameOut = [dateStamp prms.addStr '.' num2str(tetNum)];
        data = int8(zeros(216,size(wavFin,1)));
        %make header - edit 2 fields
        durInd = strcmp(header(:,1),'duration');
        header{durInd,2} = num2str(duration);
        spkInd = strcmp(header(:,1),'num_spikes');
        header{spkInd,2} = num2str(size(wavFin,1));
        %grab time base
        tBaseInd = strcmp(header(:,1),'timebase');
        tBase = sscanf(header{tBaseInd,2},'%d');
        %get wf's into right shape
        wavFin = reshape(wavFin,size(wavFin,1),200);
        
        %retransform time stamp into 4 byte
        tStamp = int32(TSFin(:,1) * tBase);
        %needs to be big endian format
        [~,~,endian] = computer;
        if strcmp(endian,'L')
            tStamp = swapbytes(tStamp);
        end
        tStamp = reshape( typecast(tStamp,'int8'),4,[]);
        
        %assign into output
        data([5:54, 59:108, 113:162, 167:216],:) = wavFin';
        data([1:4, 55:58, 109:112, 163:166],:) = repmat(tStamp,[4,1]);
        %first write new header
        fid = fopen([prms.writeDir fNameOut],'w');
        for i = 1:length(header)
            fprintf(fid,'%s %s\r\n',header{i,1},header{i,2});
        end
        fprintf(fid,'%s','data_start'); %add data marker
        fwrite(fid,data(:),'int8',0,'ieee-be'); %write data
        fprintf(fid,'\r\n%s\r\n','data_end'); %add data end marker
        fclose(fid);
        
    otherwise
        
        if prms.addNoise
            voltTrace = randi([-prms.noiseMax prms.noiseMax],prms.sRate * duration,4);
        else
            % baseline - can be all zeros
            voltTrace = int16(zeros(prms.sRate * duration,4));
        end
        
        %grab waveforms and reshape into spike sample (nSpikes*50) x tet
        %channel vector (i.e. 4)
        tempWF = permute(int16(wavFin),[2 1 3]);
        tempWF = reshape(tempWF,size(wavFin,1)*50,[]);
        %make up time indices for spikes - assuming standard AXONA format - 1ms
        %sample/spike with 200micros pre and 800 micros post threshold
        timeInd = repmat(TSFin * prms.sRate,1,1,50);
        timeInd = reshape( permute(timeInd,[3 1 2]),size(wavFin,1)*50,[]);
        %sampling points relative to threshold crossed
        samplePoints = linspace(-0.0002 * prms.sRate,0.0008 * prms.sRate,50);
        samplePoints = repmat(samplePoints',size(wavFin,1),4); %get into 'tempWF' format
        
        %get index for each spike
        spikeTimesInd = round(timeInd + samplePoints);
        %convert to subs
        colInd = repmat(1:4,length(spikeTimesInd),1);
        colInd = colInd(:);
        rowInd = spikeTimesInd(:);
        spikeTimeSub = sub2ind(size(voltTrace),rowInd,colInd);
        %add voltages to trace
        voltTrace(spikeTimeSub) = prms.wfSign * tempWF; %add waveform voltages
        voltTrace = voltTrace';
        
        switch mode
            case 'mountainsort'
                % make index of peak positions - same format as firings output from
                % sorting call
                [maxVals,chInd] = max(double(wavFin),[],3);
                [~,peakInd] = max(maxVals,[],2);
                idx = sub2ind(size(chInd),1:size(chInd,1),peakInd');
                maxCh = chInd(idx);
                
                peakInd = (peakInd + spikeTimesInd(1:50:end,1))-1;
                
                firingsOut = [maxCh; peakInd'];
                
                %         writemda(event_times,sprintf('event_times.nt%02d.mda',tetNum),'int32');
                fNameOut = [dateStamp '_' sprintf('raw.%i.mda',tetNum)];
                writemda(voltTrace,[prms.writeDir fNameOut],'int16');
                save([prms.writeDir 'firingsOutAll.mat'],'firingsOut');
                
            case 'kilosort'
                %write file
                fNameOut = [dateStamp '_Tet_' num2str(tetNum) '.dat'];
                fid = fopen([prms.writeDir fNameOut],'w');
                fwrite(fid,voltTrace,'int16');
                fclose(fid);
            case 'jrclust'
                %write file
                fNameOut = [dateStamp '_Tet_' num2str(tetNum) '.bin'];
                fid = fopen([prms.writeDir fNameOut],'w');
                fwrite(fid,voltTrace,'int16');
                fclose(fid);
        end
        
end
dirOut = prms.writeDir;
end

