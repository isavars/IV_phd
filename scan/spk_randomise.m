
function [data]=spk_randomise(data,randType,varargin)
% Randomise spike times.
%
%       data=spk_randomise(data,'random');
%       data=spk_randomise(data,'shuffle',chunk,wrapMode);
%       data=spk_randomise(data,'wrap',minShift);
%       data=spk_randomise(data,'fixedWrap',shiftValue);
%
% 'random' mode re-assigns spikes to random times in the trial. N Spike is
% preserved.
%
% 'shuffle' uses a process designed to preserve temporal structure of spike firing:
% spikes are blocked into long chunks (e.g. 10-30s), and these chunks are
% then randomly re-arranged within the trial - times within the chunks remain 
% unchanged. Additionally (if wrapMode=1) spike times within each chunk are 
% wrapped around a random point the chunk.
%
% 'wrap' wraps all spikes around a random point in the trial. Random point
% must be minShift seconds away from trial start or end (default 10s).
%
% 'fixedWrap' wraps all spikes by a set amount (such that the ensemble spiking remains
% undisturbed, temporally, but relation to position is broken). Randomisation must
% be handled to calling function.

if strcmp(randType,'shuffle')
    chunk=varargin{1};
    wrapMode=varargin{2};
elseif any(strcmp(randType,{'wrap','fixedWrap'}))
    shiftArg=varargin{1}; 
end

% 'fixedWrap' args can either be one value, or one per trial. If one value, expand to be a 1:nTrial vector %
if strcmp(randType,'fixedWrap') && numel(shiftArg)==1
    shiftArg = repmat( shiftArg, 1, length(data.trials) );
end

for ii=1:length(data.trials)
    for jj=1:length(data.trials(ii).cells)
        switch randType
            case 'shuffle'
                randST=spikeShuffle(data.trials(ii).cells(jj).st,data.trials(ii).dur,chunk,wrapMode);
            case 'wrap'
                shift= (rand(1) * (data.trials(ii).dur-(shiftArg*2))) + shiftArg;
                randST=data.trials(ii).cells(jj).st + shift;
                randST(randST>data.trials(ii).dur) = randST(randST>data.trials(ii).dur) - data.trials(ii).dur; 
            case 'fixedWrap'
                shift = shiftArg(ii);
                randST=data.trials(ii).cells(jj).st + shift;
                randST(randST>data.trials(ii).dur) = randST(randST>data.trials(ii).dur) - data.trials(ii).dur; 
            case 'random'
                randST=rand(size(data.trials(ii).cells(jj).st)) .* data.trials(ii).dur;
                randST=sort((round(randST.*10000))./10000); % Temporal resolution to 0.1ms
        end             
        % Assign new spike times to data %
        data.trials(ii).cells(jj).st=sort( randST(~isnan(randST)) );
        % Check that no new spike times go over end of trial (potential overhang is
        % ~10ms), so don't worry about moving other spikes, just cut off any overhangers.
        data.trials(ii).cells(jj).st=data.trials(ii).cells(jj).st( data.trials(ii).cells(jj).st < length(data.trials(ii).x)*(1/data.trials(ii).sample_rate) );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [randST]=spikeShuffle(ST,trialDur,chunk,wrapMode)
%%% This sub-function performs the shuffling of the spike train times. %%%
chunkStarts=0:chunk:trialDur;
chunkLengths=diff([chunkStarts trialDur]);
if chunkLengths(end)==0
    chunkStarts=chunkStarts(1:end-1);
    chunkLengths=chunkLengths(1:end-1);
end
% Split spike times by chunk - get times relative to chunk start %%
chunkedST=repmat(NaN,length(ST),length(chunkStarts)); % Pre-allocate (spikeTime,chunk) matrix
for kk=1:length(chunkStarts)
    chunkTemp=ST(   ST>=chunkStarts(kk)  &   ST<(chunkStarts(kk)+chunkLengths(kk))  );
    chunkedST(1:length(chunkTemp),kk)=chunkTemp;
end
timesInChunk=chunkedST - repmat(chunkStarts,size(chunkedST,1),1);
% Wrap times within each chunk %%
if wrapMode
    wrapPoints=rand(size(chunkLengths)).*chunkLengths;
    wrappedST=timesInChunk - repmat(wrapPoints,size(chunkedST,1),1);
    indTemp=wrappedST<0;
    chunkLengthsXSpike=repmat(chunkLengths,size(chunkedST,1),1);
    wrappedST(indTemp)=wrappedST(indTemp) + chunkLengthsXSpike(indTemp);
    timesInChunk=wrappedST;
end
% Randomise chunk starts, and add to chunked spike times %%
randomInd=randperm(length(chunkStarts));
randChunkStarts=chunkStarts(randomInd);
% Need to deal with the shorter end chunk - all chunk starts after this
% chunk in the new arrangement need to have the appropriate amount subtracted.
indTemp=randChunkStarts>randChunkStarts(end);
randChunkStarts(indTemp)=randChunkStarts(indTemp) - (chunk-chunkLengths(end));
randST=timesInChunk + repmat(randChunkStarts,size(chunkedST,1),1);        
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     chunkStarts=0:chunk:data.trials(ii).dur;
%     chunkLengths=diff([chunkStarts data.trials(ii).dur]);
%     if chunkLengths(end)==0
%         chunkStarts=chunkStarts(1:end-1);
%         chunkLengths=chunkLengths(1:end-1);
%     end

%         %% Split spike times by chunk - get times relative to chunk start %%
%         chunkedST=repmat(NaN,length(data.trials(ii).cells(jj).st),length(chunkStarts)); % Pre-allocate (spikeTime,chunk) matrix
%         for kk=1:length(chunkStarts)
%             chunkTemp=data.trials(ii).cells(jj).st(  ...
%                 data.trials(ii).cells(jj).st>=chunkStarts(kk) & data.trials(ii).cells(jj).st<(chunkStarts(kk)+chunkLengths(kk))  );
%             chunkedST(1:length(chunkTemp),kk)=chunkTemp;
%         end
%         timesInChunk=chunkedST - repmat(chunkStarts,size(chunkedST,1),1);
%         
%         %% Wrap times within each chunk %%
%         if wrapMode
%             wrapPoints=rand(size(chunkLengths)).*chunkLengths;
%             wrappedST=timesInChunk - repmat(wrapPoints,size(chunkedST,1),1);
%             indTemp=wrappedST<0;
%             chunkLengthsXSpike=repmat(chunkLengths,size(chunkedST,1),1);
%             wrappedST(indTemp)=wrappedST(indTemp) + chunkLengthsXSpike(indTemp);
%             timesInChunk=wrappedST;
%         end
%         
%         %% Randomise chunk starts, and add to chunked spike times %%
%         randomInd=randperm(length(chunkStarts));
%         randChunkStarts=chunkStarts(randomInd);
%         % Need to deal with the shorter end chunk - all chunk starts after this
%         % chunk in the new arrangement need to have the appropriate amount subtracted.
%         indTemp=randChunkStarts>randChunkStarts(end);
%         randChunkStarts(indTemp)=randChunkStarts(indTemp) - (chunk-chunkLengths(end));
%         randST=timesInChunk +
%         repmat(randChunkStarts,size(chunkedST,1),1);