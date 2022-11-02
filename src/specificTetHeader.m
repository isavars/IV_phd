function [tetFileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,tetrode_number)
tetrodeHeader = genericTetHeader;

% Insert 'num_spikes' for given tetrode
numberSpikes = num2str(num_spikes(tetrode_number));
spikesInd = find(strcmp(tetrodeHeader(:,1),'num_spikes'));
tetrodeHeader{spikesInd,2} = numberSpikes;

% % Change tetrode number in file name to correspond to DACQ
% if strcmp('Right',trialInfo.plug)
%     tetrode_number = tetrode_number + 4; % DACQ tetrodes 5 -8
% else 
%     % do nothing, DACQ tetrodes 1 - 4
% end

%% Write to file 
tetFileName = [char(trialInfo.finalTrialName) '.' num2str(tetrode_number)];
fid = fopen([trialInfo.writeDir '\' tetFileName],'w');
for j = 1:length(tetrodeHeader)
%     if j == length(tetrodeHeader)
%         fprintf(fid,'%s %s',tetrodeHeader{j,1},tetrodeHeader{j,2});
%     else
        fprintf(fid,'%s %s\r\n',tetrodeHeader{j,1},tetrodeHeader{j,2});
%     end
end
fclose(fid);

