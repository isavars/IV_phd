function [rtn] = gss()

% GSS - Get Scan Structure. 
%
%       scanData = gss;
%
% Returns SCAn application data structure:
%
% (These are the important ones for addressing data)
%
% scanData.selData = {[]} - Cell array of strings, dataset selected in base workspace.
% scanData.selTrial = []  - Index of selected trials in selected dataset.
% scanData.selCell = []   - Index of selected cells in selected dataset.
% scanData.selMaps = []   - Selected maps. Maps Parameters struct - use
%                           RATES_MAPMATCH to recover index for each dataset.
%
% (Also these, used for running the GUI)
%
% scanData.infoDataShow = []
% scanData.infoDataSort = []
% scanData.infoTrialShow = []
% scanData.infoTrialSort = []
% scanData.infoCellShow = []
% scanData.infoCellSort = []
% scanData.settings.lastDirLoad
% scanData.settings.lastDirOpen
% scanData.settings.lastDirSave
% scanData.settings.lastDirEps

set(0,'showhiddenhandles', 'on');
hAllFigs = get(0, 'children');
hFig = findobj(hAllFigs, 'flat', 'tag', 'scanBaseWindow');
set(0,'showhiddenhandles', 'off');
if isempty(hFig);
    rtn = [];
else
    rtn = guidata(hFig);
end