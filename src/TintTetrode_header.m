function [genericTetHeader] = TintTetrode_header (duration) %TintTetrode_header(trialInfo,DACQinfo)
%% Amend tetrode header
% TETfile_header = load('TetFileHeader.mat');
% genericTetHeader = TETfile_header.tet_header;
% duration = num2str(trialInfo.trial_duration*60); %change to duration coming from dat2volts
% duration = num2str(duration); %added this 

genericTetHeader = {'trial_date','Sunday, 20 Feb 2022'; ...   %change to empty from trialInfo.hdrDate
'trial_time','18:48:22'; ... %change to empty from trialInfo.trial_time
'experimenter','TOD'; ...
'comments','Logger trial';...
'duration',duration;...       
'sw_version','1.3.4.13'; ... % hc instead of DACQinfo.sw_version
'num_chans','4'; ...
'timebase','96000 hz'; ...
'bytes_per_timestamp','4'; ...
'samples_per_spike','50'; ...
'sample_rate','48000 hz'; ...
'bytes_per_sample','1'; ...
'spike_format','t,ch1,t,ch2,t,ch3,t,ch4'; ...
'num_spikes',''}; 


% % Insert date
% dateInd = find(strcmp(genericTetHeader(:,1),'trial_date'));
% genericTetHeader{dateInd,2} = trial_date; % Add date info to header 
% 
% % Insert time
% timeInd = find(strcmp(genericTetHeader(:,1),'trial_time')); % index of trial_time in header fields
% genericTetHeader{timeInd,2} = trial_time; % Add time info to header
% 
% % Insert duration
% durInd = find(strcmp(genericTetHeader(:,1),'duration')); % index of trial_time in header fields
% genericTetHeader{durInd,2} = duration; % Add time info to header

%% Make column vector
% genericTetHeader = reshape(genericTetHeader.',[],1);
% genericTetHeader = genericTetHeader(1:end-1); % remove empty cell at end




