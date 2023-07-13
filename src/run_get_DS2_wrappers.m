function run_get_DS2_info_wrapper(sws_table, writeDir)
    %run get_DS2_info_wrapper systematically over a chosen set of sleep
    %trials in sleepData table. 

    % TO DO needs to find probe data from somewhere - might need to use a
    % spreadsheet even though all the other values are in the sleepData
    % table. h.c. for now 
    
    %load sleep data table 
    load (sws_table, 'sleepData')
    % adjustable parameters 
    num_spikes = 5; %first 5 for now but needs to be all the spikes once detection improves - can be added as a variable in sleepData 
    probe_type_1 = {'r1099'};
    probe_type_2 = {'r1133', 'r1146', 'r1151', 'r1152', 'r1166', 'r1175', 'r1176', 'r1196', 'r1207', 'r1208', 'r1241','r1242', 'r1311'};
    probe_type_3 = {'r1100', 'r1142', 'r1143', 'r1144'};

    for ii = 1:height(sleepData)
        %get path for eeg data 
        eeg_data = char(strcat(sleepData.trialpath(ii), sleepData.trialname(ii), '.eeg'));%changed here so its a character string not sure why I needed to 
        %select the corect probe type and run the wrapper 
        if any(strcmp(sleepData.ratID(ii), probe_type_1))
            probe_type = 1;
            get_DS2_info(num_spikes, eeg_data, probe_type, writeDir, sws_table);
        elseif any(strcmp(sleepData.ratID(ii), probe_type_2))
            probe_type = 2; 
            get_DS2_info(num_spikes, eeg_data, probe_type, writeDir, sws_table);
        elseif any(strcmp(sleepData.ratID(ii), probe_type_3))
            probe_type =3; 
            get_DS2_info(num_spikes, eeg_data, probe_type, writeDir, sws_table);
        end 
    end 
end 
