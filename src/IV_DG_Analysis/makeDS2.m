function makeDS2()
eeg_output = [];
for it_eegs = 1:32
    %eeg 1 file doesnt have the 1 so figure that out
    if it_eegs == 1
        eeg_struct  = load_eeg(strcat('211129f_sleepHP_LP.eeg'));
    else
        eeg_struct = load_eeg(strcat('211129f_sleepHP_LP.eeg',num2str(it_eegs)));
    end
    eeg_output = [eeg_struct.eeg;eeg_output];
end

eeg_output = transpose(reshape(eeg_output,[300000,32]));

%convert to mV 

eeg_output = eeg_output.*(1.5 ./ 128);


% plot(eeg_output(1,[500:1000]))

    all_thresh_cross = [];
    thresh = 1.2;

    absolute_eeg_output = abs(eeg_output);
    thresh_cross = absolute_eeg_output > thresh;

% for it_eeg = 1:300000
%     if (eeg_output(1,it_eeg) > thresh && eeg_output(1,it_eeg + 5) > thresh && eeg_output(1,it_eeg + 10) > thresh) || (eeg_output(1,it_eeg) < -thresh && eeg_output(1,it_eeg + 5) < -thresh && eeg_output(1,it_eeg + 10) < -thresh)
%         thresh_cross = eeg_output(1,[it_eeg-14:it_eeg+17]);
%         all_thresh_cross = [all_thresh_cross;thresh_cross];
%     end
%     
% end 


% 
% all_thresh_cross = all_thresh_cross(1:20:end,:);
% whos
% 
% for it_th = 1:length(all_thresh_cross)
%     figure()
%     plot (all_thresh_cross(it_th,:))
% 
% end