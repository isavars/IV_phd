%make tetrode files from converted dacq raw output 

fileName = dir('**/*.dat');
name = fileName.name;
name = extractBefore(name, ".");
trialInfo.finalTrialName = name;
trialInfo.writeDir = fileName.folder;
    
    %% Generate Tint .tet header
    [genericTetHeader] = TintTetrode_header(trial_duration);%replaced trialInfo,DACQinfo with duration
    [tet1FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,1);
    [tet2FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,2);
    [tet3FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,3);
    [tet4FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,4);
    [tet5FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,5);
    [tet6FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,6);
    [tet7FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,7);
    [tet8FileName] = specificTetHeader(genericTetHeader,trialInfo,num_spikes,8);
    %% Append header and create tetrode file
    fprintf('Generating tetrode files....\n')
    make_tint_file(tetrode1,tet1FileName,1,trialInfo); % removed ,DACQinfo
    make_tint_file(tetrode2,tet2FileName,2,trialInfo);
    make_tint_file(tetrode3,tet3FileName,3,trialInfo);
    make_tint_file(tetrode4,tet4FileName,4,trialInfo);
    make_tint_file(tetrode5,tet5FileName,5,trialInfo);
    make_tint_file(tetrode6,tet6FileName,6,trialInfo);
    make_tint_file(tetrode7,tet7FileName,7,trialInfo);
    make_tint_file(tetrode8,tet8FileName,8,trialInfo);
%     make_tint_file(tetrode9,tet1FileName,1,trialInfo); % removed ,DACQinfo
%     make_tint_file(tetrode10,tet2FileName,2,trialInfo);
%     make_tint_file(tetrode11,tet3FileName,3,trialInfo);
%     make_tint_file(tetrode12,tet4FileName,4,trialInfo);