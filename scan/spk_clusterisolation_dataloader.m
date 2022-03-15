function [] = spk_clusterisolation_dataloader()
% Load in Cell 0, then calculate cluster isolation, for all data.
%
%     [] = spk_clusterisolation_dataloader()
%
% If you are working on the same computer as that on which you originally loaded the data into SCAN, and no directory
% names have changed, then you can simply run this function (without arguments)
%
% If not, then you need to change the matlab current directotry such that you are in the folder directly above the
% rat folders, then run the function.
%
% LR and isoD will be added as metadata tags for each cell to the SCAn structure.

SD=gss;
for ii=1:length(SD.selData)
    data = evalin('base', SD.selData{ii});
    fprintf(1,'%s ', SD.selData{ii});
    % Check to see if this is EEG only data - if so, skip %
    if isempty(data.trials(1).cells);   continue;   end
    % Set LS to load Cell 0 %
    LS=data.load_selection;
    LS.loadSpikes = 1;
    LS.wfMode = 2;
    LS.loadPos = 0;
    LS.loadEEG1 = 0;
    LS.loadEEGAll = 0;
    for jj=1:8;
        LS.tets(jj).load0 = 1;
        LS.tets(jj).trimNull = 0;
        % Get list of which cells were loaded - if use trimNull, not all 0 spikes necessarily go to cell 0.
        tetInd = find(cat(1,data.trials(1).cells.tet)==jj);
        LS.tets(jj).cellList = cat(1,data.trials(1).cells(tetInd).cellnum);
        if ~isempty(tetInd)
            cutEnd=fliplr(strtok(fliplr( data.trials(1).cells(tetInd(1)).cut ),'_' )); % cut after last '_'
            LS.tets(jj).cutTag=['_' strtok(cutEnd,'.')];
            LS.tets(jj).cut = '';
        end
    end
    % Check if path to DACQ files specified in data exists .. %
    currDir=pwd;
    try cd(LS.pathName)
        cd(currDir);  % If it does, do nothing here, pass path to LOAD_MAIN, below. 
    catch
        ratDirName = fliplr(  strtok(  fliplr(LS.pathName),  '\/'  )   );   % Otherwise, assume that the matlab currect directory is the directory
        LS.pathName = [currDir '\' ratDirName];                 % that contains the rat directories.
    end
    % Load in data with cell 0 %    
    dataC0 = load_main([], LS, LS.trialNames, LS.pathName);
    % Run Cluster Isolation analysis %
    for jj=1:length(data.trials)
        tetList=unique(cat(1,dataC0.trials(jj).cells.tet));     % tetList: which tetrodes are in use
        for kk=1:length(tetList)
            tetInd = find(   cat(1, dataC0.trials(jj).cells.tet)  ==  tetList(kk)   );  % tetInd: numerical index into data.trials.cells, which cells are from the tet in question
            cellList = cat(1, dataC0.trials(jj).cells( tetInd ).cellnum);               % cellList: list of cells (cellnum, i.e. TINT numbers) on the tetrode in question.
            % Convert the amp values in the SCAn structure into a single array %
            amps = [];   cellIndForSpikes = [];
            for mm=1:length(tetInd)
                if ~isempty(dataC0.trials(jj).cells( tetInd(mm) ).st)
                    amps = cat(1, amps, dataC0.trials(jj).cells( tetInd(mm) ).wf_amps);
                    cellIndForSpikes = cat(1, cellIndForSpikes, ones(length(  dataC0.trials(jj).cells( tetInd(mm) ).st  ),1).*cellList(mm));
                else
                    % If no spikes, put in a single dummy spike, this will be automatically set to NaN, as Mahanalobis distance only defined if n spikes greater than d.f.
                    amps = cat(1, amps, [0 0 0 0]);
                    cellIndForSpikes = cat(1, cellIndForSpikes, cellList(mm));
                end
            end
            [LR, IsoD] = spk_clusterisolation(amps, cellIndForSpikes);
            % Assign cluster data to normal data cells.user %
            tetInd = find(cat(1, data.trials(jj).cells.tet)==tetList(kk));
            for mm=1:length(tetInd)
                data.trials(jj).cells(tetInd(mm)).user(1).LR = LR(mm);
                data.trials(jj).cells(tetInd(mm)).user(1).IsoD = IsoD(mm);
            end
        end
    end
    assignin('base', SD.selData{ii}, data);
    clear data
end
