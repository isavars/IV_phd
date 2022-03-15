function [] = analysis_clusterisolation(SD)
% Load in Cell 0, then calculate cluster isolation, for all data.

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
        tetInd = find(  cat(1,data.trials(1).cells.tet)  ==jj);
        LS.tets(jj).cellList = cat(1,data.trials(1).cells(tetInd).cellnum);  % Get list of which cells were loaded (if some cells from cut not loaded should go into '0' cluster) %
        if ~isempty(tetInd)
            cutEnd=fliplr(strtok(fliplr( data.trials(1).cells(tetInd(1)).cut ),'_' )); % cut after last '_'
            LS.tets(jj).cutTag=['_' strtok(cutEnd,'.')];    % For loading multiple trials/tets, need the cut Tag, not the whole cut name.
            LS.tets(jj).cut = '';
        end
    end
    
%     if strcmp(LS.pathName(end),'\')
%         pathName = [rootDir LS.pathName(end-5:end)];
%         data.load_selection.pathName = data.load_selection.pathName(1:end-1);
%     else
%         pathName = [rootDir LS.pathName(end-4:end)];
%     end

    dataC0 = load_main([], LS, LS.trialNames, LS.pathName);
    
    %%% Run Cluster Isolation analysis %%%
    for jj=1:length(data.trials)
        tetList=unique(cat(1,dataC0.trials(jj).cells.tet));
        for kk=1:length(tetList)    % For each tetrode
            currTetInd = cat(1, dataC0.trials(jj).cells.tet)==tetList(kk);  % Get list of cells on current tet
            currCells=dataC0.trials(jj).cells(currTetInd);
            amps = [];     cellNums = [];
            for mm=1:length(currCells)
                amps = [amps; double(currCells(ii).wf_amps)];
                cellNums = [cellNums; repmat(currCells(ii).cellnum, length(currCells(ii).st), 1)];
            end
            [LR, IsoD] = spk_clusterisolation(amps,cellNums);
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
    pack
end

helpdlg('Cluster Isolation values have been successfully calculated. They are shown in the ''Cells''', ...
         'properties window. To output to excel, use ''Analyse'' -> ''Output 2 Excel''' ,'Cluster Isolation Done');

