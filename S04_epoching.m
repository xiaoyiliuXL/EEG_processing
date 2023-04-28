%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------
% Goal of the script :
% Epoch the time series, and reject bad epoch
% Sliding window method (in erplab) is used to detect bad epoches
% 
% isVisualInspect = 1 (recommended, to do a visual inspection)
% When visually inspecting the epoches, you will be prompted:
% 1. Don't reject: (epoch ids that you think is a good epoch, 
%                   but the algorithm thinks is a bad epoch)
% 2. Reject: (epoch ids that you think is a bad epoch, 
%             but the algorithm thinks is a good epoch)
%
% ----------------------------------------------------------------------
% Template created by Xiaoyi LIU (xiaoyi.x.liu@gmail.com)
% Last update : April 2023
% Project : -
% Version : -
% ----------------------------------------------------------------------

%% -------------------- PREPARATION -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% add eeglab toolbox to path
addpath('XXXX'); 
% e.g., '/Users/xl4251/Documents/MATLAB/eeglab2021.1'

% init eeglab (and functions in the toolbox)
eeglab;
rng(6); % a random seed

%% -------------------- INITIATING -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing parameters to specify
subjNums = [1]; % id for all VALID participants
nSubj = length(subjNums); % the total number of VALID participantas

% a file that has the event triggers for each condition you are intereseted
% in
binlisterFile = 'binlister.txt';

% epoch, making sure only trials with an end are included
tEndTrigger = {'S 12'}; % the trigger for trial-end                         
timeWin = [-5 0.5]; % time window of the epoch around the time of trial-end

% epoching, relevant to a certain event
timeWin_toEvent = [-0.5 1]; % time window around the time of the event

% detect bad epoch with moving window
threshold = 100; % amplitude
toi = [-0.5 1]; % time window of interest
winSize = 200; % window size
winStep = 50; % window sliding step

isVisualInspect = 1; % whether or not to do a visual inspection

% channels
nonscalp_channels = {};
EOG_channels = {};
EOG_channels_mat = cell2mat(EOG_channels);
scalp_channels = 1:32;

isPlot = 1; % plot or not after processing for each subject

%% -------------------- SUBJ LOOP -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSubj = 1:nSubj
    
    % input files
    curr_subjno = subjNums(iSubj); 
    curr_inputname = sprintf('S%02d_icareject.set', curr_subjno);
    
    % output files
    output_dir = fullfile('results', sprintf('S%02d', curr_subjno));
    curr_outputname = sprintf('S%02d_epoch.set', curr_subjno);
    
    if exist(fullfile(output_dir, curr_outputname), 'file')
        continue
    else

        %% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ICA_EEG = pop_loadset('filename', curr_inputname, 'filepath', output_dir);
        
        %% epoching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % epoch, making sure only trials with an end are included
        ICA_EEG = pop_epoch(ICA_EEG, tEndTrigger, timeWin, 'epochinfo', 'yes');

        % remove eye channels
        erpEPOCH_EEG = pop_select(ICA_EEG, 'nochannel', cell2mat([nonscalp_channels, EOG_channels]));

        % epoching w erplab
        erpEPOCH_EEG  = pop_creabasiceventlist(erpEPOCH_EEG, 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' } ); 
        
        erpEPOCH_EEG  = pop_binlister(erpEPOCH_EEG, 'BDF', binlisterFile, 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG' );

        erpEPOCH_EEG = pop_epochbin(erpEPOCH_EEG, timeWin_toEvent,  'pre'); %(baseline removed)

        % detect bad epoch with moving window
        erpEPOCH_EEG  = pop_artmwppth(erpEPOCH_EEG, 'Channel',  scalp_channels, 'Flag',  1, ...
                        'Threshold',  threshold, 'Twindow', toi, ...
                        'Windowsize',  winSize, 'Windowstep', winStep ); 
        
        %% visual inspection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isVisualInspect
            
            assignin('base', 'EEG', erpEPOCH_EEG);

            pop_eegplot(erpEPOCH_EEG, 1, 1, 1);

            no_reject = []; reject = [];
            no_reject_input = strsplit(input('Don''reject: ','s'));
            
            for i = 1:size(no_reject_input,2)
                no_reject(i) = str2num(no_reject_input{i});
            end

            reject_input = strsplit(input('Reject: ', 's'));
            
            for i = 1:size(reject_input,2)
                reject(i) = str2num(reject_input{i});
            end

            comp_handles = findobj('-regexp', 'tag', '^EEGPLOT.*');
            while 1
                try
                    waitforbuttonpress
                catch
                    if ~any(isgraphics(comp_handles))
                        break
                    end
                end
            end

        end
        
        erpEPOCH_EEG.reject.rejmanual = evalin('base', 'ALLEEG(end).reject.rejmanual');

        %% remove bad epoches %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        nEpoch = length(erpEPOCH_EEG.epoch);
        
        for iEpoch = 1:nEpoch
            if any(iEpoch == no_reject)
                erpEPOCH_EEG.epoch(iEpoch).eventflag = 0;
                erpEPOCH_EEG.reject.rejmanual(iEpoch) = 0;
            elseif any(iEpoch == reject)
                erpEPOCH_EEG.epoch(iEpoch).eventflag = 1;
                erpEPOCH_EEG.reject.rejmanual(iEpoch) = 1;
            end
%             badepoch(iEpoch) = erpEPOCH_EEG.epoch(iEpoch).eventflag;
        end
        
        badepoch_idx = find(erpEPOCH_EEG.reject.rejmanual);
        
        % remove bad epoches
        EPOCHRM_EEG = pop_rejepoch(erpEPOCH_EEG, erpEPOCH_EEG.reject.rejmanual);

        %% saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        pop_saveset(EPOCHRM_EEG, 'filename', fullfile(output_dir, curr_outputname));
        
        %% documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        epoch_filename = sprintf('epoch_%02d.txt', curr_subjno);
        epoch_file = fopen(fullfile(output_dir, epoch_filename), 'w');

        fprintf(epoch_file, 'Trial indices rejected after ICA: ');
        fprintf(epoch_file, '%i\t', badepoch_idx);
        fprintf(epoch_file, '\n\nNumber of trials rejected: ');
        fprintf(epoch_file, '%2d', sum(erpEPOCH_EEG.reject.rejmanual));

        fclose('all');
        
    end
end