
function [rtn] = load_tet(trialName, LS, tetN, setFile, pos)

% Read a tetrode file and a cut file together and return a structure array containing
% the time and waveform data for each TINT defined cell on that tetrode.
%
%        C = load_tet('trialName', LS, tetN, setFile);
%         
% LS is the output of LOAD_SELECTDATA, tetN is the tet number to load, 
% setFile is the output of LOAD_SET.
%
% C is a .cells struct:
%
%   .cells.cellnum = [];
%   .cells.tet = [];
%   .cells.cut = [];
%   .cells.scalemax = [];
%   .cells.user = struct([]);
%   .cells.st = [];
%   .cells.wf_all = [];
%   .cells.wf_means = [];
%   .cells.wf_amps = [];
%   .cells.wf_mode = 0;


%%% Read Spike Files %%%
spikes = load_spikes([trialName, '.', num2str(tetN)]);

%%% Read Cut file %%%
if isempty(LS.tets(tetN).cut)
    cutFileName = [trialName, LS.tets(tetN).cutTag, '.cut'];
else
    cutFileName = LS.tets(tetN).cut;
end
fid = fopen(cutFileName,'r');
if fid==-1
    % Cut file not found: depending on caller input, can either throw an error or return a structure with no spikes %
    if ~LS.skipMissingTets
        error(['SCAN Error: Cut file ', cutFileName, ' not found.']);
    else
       disp(['SCAN Warning, dataset ' LS.dataName ': Cells on cut file ', cutFileName, ' have been skipped']);
       cut = repmat(-1, size(spikes.times)); % Most convenient not to declare the empty structure directly, but to make a fake cut file, in which none of the specified cells will be found (all cells cluster ID = -1).
    end
else
    % Read cut file %
    n=0;                                                    %
    while 1                                                 %
        tline = fgetl(fid);                                 % Look for start of exact cut data.
        if strmatch('Exact_cut_for',tline), break, end      % 
        n=n+1;                                              % The final returned tline is the header for this.
    end                                                     %
    dataStart = n;                                          %
    fclose(fid);                                            
    cut = textread(cutFileName,'%n','headerlines',dataStart+1,'delimiter',' ');       % Read spike cluster data
    % Tetrode-Cut Mismatch %%    
    if length(spikes.times) ~= length(cut)
        error(['Trial ', trialName, ': Number of spikes in cut and tetrode file don''t match. Cut file saved with wrong trial?']);
    end
end
% Remove from cut file reference to cells overhanging trial end %%
cut = cut(spikes.times < size(pos.led_pos,1)/pos.sample_rate); % Use nPosSamp/SR: stops spikes overhanging the rounding error of DACQ1 SR.

%%% Assign spikes to cells %%%
if LS.tets(tetN).load0
    LS.tets(tetN).cellList = [LS.tets(tetN).cellList; 0];    % Put cell 0 at end of list first, if necessary.
end
for ii=1:length(LS.tets(tetN).cellList) % For each cell
    % Get index of spikes belonging to this cell %%
    if LS.tets(tetN).cellList(ii)==0
        % Load spike 0: Make special index, to load all non-selected cells to cell 0 %
        nsCells = unique([0, setdiff( (0:30), LS.tets(tetN).cellList )]);
        currCellInd = [];
        for jj=1:length(nsCells)
            currCellInd = [currCellInd; find(cut==nsCells(jj))];
        end
    else
        % Otherwise just look for index of cell number %
        currCellInd = find(cut == LS.tets(tetN).cellList(ii));
    end

    % Select spike data for this cell %%
    rtn(ii).st = spikes.times(currCellInd);
    rtn(ii).wf_all = [];    %
    rtn(ii).wf_means = [];  % Defaults for wf data (wfMode=0, ie no wf data)
    rtn(ii).wf_amps = [];   %
    switch LS.wfMode
        case 2 % Only amps and mean wfs. 2 for backwards compat.
            spkTemp = double(spikes.waves(currCellInd,:,:));
            if size(spkTemp,1)~=0
                rtn(ii).wf_means = single(squeeze(mean(spkTemp)));
                rtn(ii).wf_amps = uint8( squeeze(max(spkTemp,[],2)) - squeeze(min(spkTemp,[],2)) );
                if length(rtn(ii).st)==1
                    rtn(ii).wf_amps=rtn(ii).wf_amps';
                end                    
            end
        case 1 % All waveform data. 1 for backwards compat.
            rtn(ii).wf_all = spikes.waves(currCellInd,:,:);
    end
    % Other cell metadata %%
    rtn(ii).cut = cutFileName;
    rtn(ii).cellnum = LS.tets(tetN).cellList(ii);
    rtn(ii).tet = tetN;
    rtn(ii).scalemax = setFile.fullscale(tetN,:);
    rtn(ii).wf_mode = LS.wfMode;
    rtn(ii).user.brain_region = '';
end





