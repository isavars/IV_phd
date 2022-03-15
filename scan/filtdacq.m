
function [] = filtdacq()

% Filter spikes by changing threshold levels in DACQ files.

% Open figure + set callbacks %
hFig = open('filtdacq.fig');
hApp = guihandles(hFig);
% NB. Callbacks are set here (not in fig file), as callback function handles 
% need to be created in a place where callack functions are visible.
set(hApp.open, 'callback', @filtdacqOpen);
set(hApp.save, 'callback', @filtdacqSave);
set(hApp.loadCut, 'callback', @loadCutFile);
set([hApp.thr_1 hApp.thr_2 hApp.thr_3 hApp.thr_4], 'callback', @filtdacqUpdateApp);
set([hApp.dotsize_big hApp.dotsize_small], 'callback', @filtdacqDotsize);

% Cluster Colours %
clusterCMap = [0.7 0.7 0.7; 0 0 1; 0 1 0; 1 0 0; 1 1 0; 0 1 1; 0 0.5 0];

% Initialise data vars %
bindata = [];           % Complete tetrode file, in format int8 = byte
spikedata = [];         % Spike records (including timestamp), format (1:nspike, 1:54, 1:4)
spikeamps = [];         % Spike amplitudes, format (1:nspike, 1:4)
spikerms = [];          % Spike energies, format (1:nspike, 1:4)
spikemax = [];          % Spike peaks, format (1:nspike, 1:4)
scalemax = [0 0 0 0];   % In microvolts
nspike = 0;             % This refers to the original file, not the filtered dataset.
thr = [64 64 64 64];    % Threshold, as absolute int8 sample number.
filtpass = [];          % Index of spikes that are filtered IN.
cut = [];               % 1:nSpike index of cluster assignments
fileheader = '';

filename = '';      % This will include '.1' etc for the tet number
datadir = pwd;      % This will be set to loading dir after first data load
workingdir = pwd;   % This is for reference, to cd back to.

    %%% Callbacks %%%    
    function [] = filtdacqOpen(hObject, eventData)
    %%% Open a tetrode file, get amplitudes and assign to gui vars %%%
        % File stuff %
        cd(datadir);
        [filename, pathname, filterindex] = uigetfile({'*.1'; '*.2'; '*.3'; '*.4'; '*.5'; '*.6'; '*.7'; '*.8'}, 'Choose a tetrode file to load');
        fid = fopen([pathname '\' filename],'r','ieee-be');        % 'ieee-be' is machine format: needs to be big endian
        bindata = fread(fid,'int8');     % 'int8': each byte => integer -127:128.
        fclose(fid);
        fid = fopen([pathname '\' filename(1:end-1) 'set']);
        settxt = textscan(fid, '%s');
        fclose(fid);
        datadir = pathname;
        cd(workingdir);
        % Process spike file into waveforms, get amplitudes %
        fileheader = char(abs(bindata(1:400)))';              
        m = (findstr('num_spikes ',fileheader))+11; 
        nspike = str2double(fileheader(m:(m+9)));              
        ds = (findstr('data_start',fileheader))+10;        % Data start
        spikedata = bindata(ds:(ds+(nspike*216)-1));           
        spikedata = shiftdim(reshape(spikedata,54,4,nspike), 2);
        filtpass = 1:nspike;
        cut = zeros(1,nspike);
        fileheader = fileheader(1:ds-1);
        % Amps %
        spikemax = squeeze(max(spikedata(:,5:54,:), [], 2));
        spikemin = squeeze(min(spikedata(:,5:54,:), [], 2));
        spikeamps = spikemax - spikemin;
        % RMS %
        spikerms = sqrt(mean( squeeze(spikedata(:,5:54,:).^2), 2 ));
        % Get Gains + Thr from set file %
        for ii=1:4
            ind = strmatch(['gain_ch_' num2str((filterindex-1)*4 + ii-1)], settxt{1}, 'exact');
            scalemax(ii) = str2double(settxt{1}{ind+1});
        end
        ind = strmatch('ADC_fullscale_mv', settxt{1});
        scalemax = (str2double(settxt{1}{ind+1}) ./ scalemax) .* 1000;
        ind = strmatch('threshold ', settxt{1});
        thr = repmat( floor((str2double(settxt{1}{ind+1})/32767) .* 127), 1, 4);
        % Update app figure %
        set(hApp.tetrode_name, 'string', filename);
        set(hFig, 'name', ['filtdacq: ' filename]);
        set([hApp.thr_1 hApp.thr_2 hApp.thr_3 hApp.thr_4], 'value', ((thr(1)/127)-0.5)/0.5); % Take thr from ch1 only.
        filtdacqUpdateApp;            
    end
    
    function [] = loadCutFile(hObject, eventData)
    %%%%%%%%%%%%%%%%%%%%%%% Load a cut file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cd(datadir);
        [cutfilename, pathname] = uigetfile('*.cut', 'Choose a cut file');
        cd(workingdir);
        fid = fopen([pathname '\' cutfilename]);
        tempText=textscan(fid,'%s');
        fclose(fid);
        dataStart = find(strcmp('Exact_cut_for:',tempText{1}),1);
        cut = str2double( tempText{1}(dataStart+4 : end) );        
        % Tetrode-Cut Mismatch %%    
        if nspike ~= length(cut)
            error('Number of spikes in cut and tetrode file don''t match. Cut file from wrong trial?');
        end
    end
    
    function [] = filtdacqSave(hObject, eventData)
    %%% Save the tetrode file with the currently filtered set of spikes %%%
        % Recreate filtered binary spike data %
        filtspikedata = spikedata(filtpass, :, :);
        filtspikedata = reshape(filtspikedata, size(filtspikedata,1), 216);
        filtspikedata = reshape(filtspikedata', 1, size(filtspikedata,1)*216);
        % Alter the header %
        nspikenew = [num2str(length(find(filtpass))) '          '];
        m = (findstr('num_spikes ',fileheader))+11;
        newfileheader = fileheader;
        newfileheader(m:m+9) = nspikenew(1:10);
        newfileheader = int8(newfileheader);
        % Put together new tet file %
        newbindata = [newfileheader filtspikedata [13 10 100 97 116 97 95 101 110 100 13 10]];  % Tet file termination sequence
        % Write new tetrode file %
        cd(datadir);
        fid = fopen([filename(1:end-2) '_thr.' filename(end)], 'w', 'ieee-be');
        fwrite(fid, newbindata, 'int8');
        fclose(fid);
        % Copy other data files to _thr name %
        ext = {'pos', 'eeg', 'set', 'inp'};
        for ii=1:length(ext)
            fid = fopen([filename(1:end-1) ext{ii}], 'r', 'ieee-be');
            if fid==-1;  continue;  end   % If INP file doesn't exist.
            tempbindata = fread(fid);
            fclose(fid);
            fid = fopen([filename(1:end-2) '_thr.' ext{ii}], 'w', 'ieee-be');
            fwrite(fid, tempbindata);
            fclose(fid);
        end       
        cd(workingdir);
    end        

    function [] = filtdacqThrLock(hObject, eventData)
    %%% Turns the b, c + d thr sliders on or off in response to 'lock' checkbox %%%   
        if get(hApp.thr_lock, 'value')
            set([hApp.thr_2 hApp.thr_3 hApp.thr_4], 'enable', 'off');
        else
            set([hApp.thr_2 hApp.thr_3 hApp.thr_4], 'enable', 'active');
        end
    end

    function [] = filtdacqDotsize(hObject, eventData)
    %%% Controls dotsize radio group, and calls filtdacqUpdateApp to change dots %%%
        set([hApp.dotsize_big hApp.dotsize_small], 'value', 0);
        set(hObject, 'value', 1);
        filtdacqUpdateApp;
    end

    function [] = filtdacqUpdateApp(hObject, eventData)
    %%% Update the gain and thr fields, then filter the spikes and plot them %%%    
        % Set plot mode to amps or RMS %
        if get(hApp.ampSelect,'value')
            amps_mode = 1;
        elseif get(hApp.rmsSelect,'value')
            amps_mode = 0;
        end
        % If gains locked, make them all move together %
        if get(hApp.thr_lock, 'value')
            sl = get(hApp.thr_1, 'value');
            set([hApp.thr_2 hApp.thr_3 hApp.thr_4], 'value', sl);
        end
        % generate and display thresholds + gains %
        sl = get([hApp.thr_1 hApp.thr_2 hApp.thr_3 hApp.thr_4], 'value');
        thr = floor((( cat(1,sl{:}) .* 0.5 ) + 0.5)' .* 127);
        for ii=1:4
            set( hApp.(['thr_abs_' num2str(ii)]), 'string', num2str( (thr(ii)/127)*scalemax(ii), '%3.1f' ) );
            set( hApp.(['scalemax_' num2str(ii)]), 'string', num2str(round(scalemax(ii))));
        end
        % Filter Amps %
        thr_array = repmat(thr, nspike, 1);
        filt_ch = spikemax > thr_array;
        filtpass = (sum(filt_ch,2)) > 0;
        % Plot Amps %
        for ii=1:3  % Axes y
            for jj=(ii+1):4  % Axes x
                axes(hApp.(['axes_' num2str(jj) '_' num2str(ii)]));
                delete(get(gca, 'children'));
                if amps_mode
                    plot_x = spikeamps(filtpass,jj);
                    plot_y = spikeamps(filtpass,ii);
                    set(gca,'xlim',[0 255],'ylim',[0 255]);
                else
                    plot_x = spikerms(filtpass,jj);
                    plot_y = spikerms(filtpass,ii);
                    set(gca,'xlim',[10 110],'ylim',[10 110]);
                end
%                 if get(hApp.dotsize_big, 'value')
%                     plot_x = [plot_x; plot_x+1; plot_x; plot_x+1]; % Try to simulate TINT dots by
%                     plot_y = [plot_y; plot_y; plot_y+1; plot_y+1]; % making 4 pixels per dot 
%                 end
                cutPlot = cut(filtpass);
                cellList=unique(cutPlot);
                for kk=1:length(cellList)
                    cellInd=cutPlot==cellList(kk);
                    if cellList(kk)==0;  
                        m=1;    
                    else
                        m=4;    
                    end
                    line('xdata',plot_x(cellInd),'ydata',plot_y(cellInd),'linestyle','none','marker','.','markersize',m,'markeredgecolor',clusterCMap(cellList(kk)+1,:),'markerfacecolor','none');
                end
            end
        end                
    end



end
        
        
        
        
        
        
        
        
        
        
        
        
        
        