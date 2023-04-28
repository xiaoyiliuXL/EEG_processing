%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------
% Goal of the script :
% Preprocessing steps: bad channel interpolation
% 
% Bad channels are defined as channels whose variability is 
% XX (e.g., 3) sd away from the mean of all channels
%
% Spherical interpolation is used.
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

badchan_zthresh = 3;

% channels
scalp_channels = 1:32;

isPlot = 1; % plot or not after processing for each subject

%% -------------------- SUBJ LOOP -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSubj = 1:nSubj

    % input files
    curr_subjno = subjNums(iSubj); 
    curr_inputname = sprintf('S%02d_filtered.set', curr_subjno);
    curr_inputname_ica = sprintf('S%02d_filtered_ica.set', curr_subjno);
    
    % output files
    output_dir = fullfile('results', sprintf('S%02d', curr_subjno));
    curr_outputname = sprintf('S%02d_nobadchan.set', curr_subjno);
    curr_outputname_ica = sprintf('S%02d_nobadchan_ica.set', curr_subjno);

    %% load filtered data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FILTERED_EEG = pop_loadset('filename', curr_inputname, 'filepath', output_dir);
    FILTERED_ICA_EEG = pop_loadset('filename', curr_inputname_ica, 'filepath', output_dir);

    %% detect bad channels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dat = FILTERED_EEG.data(scalp_channels, :);
    dat = zscore(dat(:));
    dat = reshape(dat, [length(scalp_channels), FILTERED_EEG.trials * FILTERED_EEG.pnts]);

    channel_std = std(dat, [], 2);
    bad_channels = find(channel_std > badchan_zthresh);
    
    fprintf('Bad channels: %i\t', bad_channels);

    %% interpolate bad channels %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    INTERP_EEG = eeg_interp(FILTERED_EEG, bad_channels, 'spherical');
    INTERP_ICA_EEG = eeg_interp(FILTERED_ICA_EEG, bad_channels, 'spherical');
    
    %% saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pop_saveset(INTERP_EEG, 'filename', fullfile(output_dir, curr_outputname));
    pop_saveset(INTERP_ICA_EEG, 'filename', fullfile(output_dir, curr_outputname_ica));
    
    %% documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    bad_chan_filename = sprintf('bad_chan_%02d.txt', curr_subjno);
    bad_chan_file = fopen(fullfile(output_dir, bad_chan_filename), 'w');
    fprintf(bad_chan_file, 'Identified bad channels: ');
    fprintf(bad_chan_file, '%i\t', bad_channels);

    fclose('all');
    
    %% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isPlot
        
        pop_eegplot(FILTERED_EEG, 1, 1, 1); %EEG
        pop_eegplot(INTERP_EEG, 1, 1, 1); %EEG

        comp_handles = findobj('-regexp', 'tag', '^eegplot.*');
        while 1 %
            try
                waitforbuttonpress
            catch
                if ~any(isgraphics(comp_handles))
                    break
                end
            end
        end % will end when the graph is closed
        
    end
    
end
