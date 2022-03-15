function [results] = analysis_output2excel(SD)
% Ouput analysis results to Excel.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Call GUI to get user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig = figure('units', 'pixels', 'position', [(SD.screenPix(3)-500)/1.3 (SD.screenPix(4)-500)/1.3 500 500], ...
               'menubar','none','numbertitle','off','name','Select analyses to run ..','closerequestfcn','delete(gcf)','renderer','painters');         
hPanAnalysis = uipanel('units','pixels','position',[5 5 240 480],'title','Analyses ..');
hPanFilter = uipanel('units','pixels','position',[255 200 240 250],'title','Filter by Trial ..');
uicontrol('units', 'pixels', 'position', [275 100 120 35], 'string', 'OK', 'callback', 'uiresume');
% Filter by trial %
filtLetterPos=[0.1 0.65 0.8 0.15];   
uicontrol('parent',hPanFilter,'units','normalized','position',filtLetterPos+[0 0.08 0 0],'style','text','string','Filter by last letter of trial name ..');
uicontrol('parent',hPanFilter,'units','normalized','position',filtLetterPos,'style','edit','tag','trialFilterByLetter','BackgroundColor','w');
filtNumberPos=filtLetterPos-[0 0.3 0 0];
uicontrol('parent',hPanFilter,'units','normalized','position',filtNumberPos+[0 0.08 0 0],'style','text','string','Filter by trial number ..');
uicontrol('parent',hPanFilter,'units','normalized','position',filtNumberPos,'style','edit','tag','trialFilterByNumber','BackgroundColor','w');
uicontrol('parent',hPanFilter,'units','normalized','position',filtNumberPos+[0 -0.3 0 0],'style','text','string','(Use commas for multiple trials e.g. a,b,c or 1,2,3)');
%%%                               %%%
%%% This defines list of analyses %%%
%%%                               %%%
              % 'Analysis Name','Function Name'; ...
analysisList = {'Theta Modulation','thetaMod'; ...
                'Complex Spike Index','complexSpikeIndex'; ...
                'Cluster isolation','clusterIsolation'};
            
yStep=0.1;  yPos=0.1;
for ii=1:length(analysisList)
    uicontrol('parent',hPanAnalysis,'units','normalized','position',[0.1 1-yPos 0.8 0.04],'value',0,'style', 'checkbox', 'string', analysisList{ii,1},'horizontalalignment','left','tag',analysisList{ii,2});
    yPos=yStep+yPos;
end
%----%
uiwait
%----%
hUI = guihandles(hFig);
userInput.analToRun={};
userInput.analysisLabels={};
for ii=1:length(analysisList)
    if get(hUI.(analysisList{ii,2}), 'value');
        userInput.analToRun{end+1}=analysisList{ii,2};
        userInput.analysisLabels{end+1}=analysisList{ii,1};
    end
end
userInput.trialFilterByLetter=get(hUI.trialFilterByLetter,'string');
userInput.trialFilterByNumber=str2num(get(hUI.trialFilterByNumber,'string')); %#ok<ST2NM>
close(hFig);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Loop through the selected data, calling analysis sub-functions, then export to Excel %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count=1;
for ii=1:length(SD.selData)
    data = evalin('base', SD.selData{ii});
    disp(SD.selData{ii});
    for jj=1:length(data.trials)
        % Filter for trial environment %
        if ~isempty(userInput.trialFilterByNumber)
            if ~any( userInput.trialFilterByNumber == jj )
                continue
            end
        elseif ~isempty(userInput.trialFilterByLetter)
            if isempty( strfind(userInput.trialFilterByLetter, data.trials(jj).trialname(end) ))
                continue
            end
        end
            
        %%% Get properties by cell %%%
        for kk=1:length(data.trials(jj).cells)
            %%% Feed all relevant data to the required sub-function %%%
            for mm=1:length(userInput.analToRun)
                dataStruct.data=data;
                dataStruct.trialIndex=jj;
                dataStruct.cellIndex=kk;
                fHandle=str2func(userInput.analToRun{mm});
                results.(userInput.analToRun{mm})(count) = fHandle(dataStruct);
            end
            % Metadata %
            results.dataset{count} = SD.selData{ii};
            results.trial{count} = data.trials(jj).trialname;
            results.trial_number(count) = jj;
            results.tet(count) = data.trials(jj).cells(kk).tet;
            results.cell(count) = data.trials(jj).cells(kk).cellnum;
            count = count + 1;
        end
    end
end
%%%% Export to Excel %%%%
% First convert the structure of different results into a N_Results x N_Measures cell string array 
resultsArray = cellstr(num2str(results.(userInput.analToRun{1})'));
for ii=2:length(userInput.analToRun)
    resultsArray = [resultsArray cellstr(num2str(results.(userInput.analToRun{ii})'))];
end
% Add metadata at the left-hand side %
resultsArray = [results.dataset', results.trial', cellstr(num2str(results.trial_number')), cellstr(num2str(results.tet')), cellstr(num2str(results.cell')), resultsArray];
% Add column labels at the top %
resultsArray = [{'Dataset','Trial','Trial Number', 'Tet', 'Cell', userInput.analysisLabels{1:end}}; resultsArray];
% Write array to excel sheet %
xlswrite('SCAn Output', resultsArray, 1, 'A1');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Below are all the actual data analysis function calls (each one within a sub-function) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rtn]=thetaMod(dataStruct)
% Theta Modulation %
[~, ~, tm, ~]=spk_intrfreqautocorr(dataStruct.data.trials(dataStruct.trialIndex),'analyseEEG',0,'cellIndex',dataStruct.cellIndex);
rtn=tm.cells(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rtn]=complexSpikeIndex(dataStruct)
% Complex Spike Index %
[~, rtn]=spk_csi(dataStruct.data.trials(dataStruct.trialIndex).cells(dataStruct.cellIndex));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rtn]=clusterIsolation(dataStruct)
% Cluster Isolation %
rtn=dataStruct.trialIndex;




