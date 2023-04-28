%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------
% Goal of the script :
% Preprocessing steps: re-reference, downsampling, hi and lo filter
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

% define the filepath of the channel map
chanmap_path = 'XXXXX';
% e.g., '/Users/xl4251/Documents/MATLAB/eeglab2021.1/plugins/dipfit/standard_BEM/elec/standard_1005.elc';

% init eeglab (and functions in the toolbox)
eeglab;

%% -------------------- INITIATING -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing parameters to specify
subjNums = [1]; % id for all VALID participants
nSubj = length(subjNums); % the total number of VALID participantas

ref = []; %     reference:      []           =       convert to average reference
          %                     [int vector] =       new reference electrode number(s)
          %                     'Cz'         =       string
          %                     {'P09' 'P10} =       cell array of strings

srate = 250; % sampling rate (downsampling)

locutoff = 0.1; 
hicutoff = 50;
locutoff_ica = 1;
hicutoff_ica = 50;

isPlot = 1; % plot or not after processing for each subject

% Define current filepath
curr_path = fileparts(which("S01_filtering.m"));
rawdata_dir = fullfile(curr_path, '/EEG_data/');

%% -------------------- SUBJ LOOP -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSubj = 1:nSubj

    % input file
    curr_subjno = subjNums(iSubj); 
    curr_setname = ['S', curr_subjno, '_EP.vdhr']; % should be the way you code the filename
    
    % output files
    output_dir = fullfile('results', sprintf('S%02d', curr_subjno));
    output_name = sprintf('S%02d_filtered.set', curr_subjno); % filtered data for analysis
    output_name_ica = sprintf('S%02d_filtered_ica.set', curr_subjno); % filtered data for ICA

    if ~exist(output_dir, 'dir') % create subject output folder if haven't
        mkdir(output_dir);
    end

    %% load raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % raw data
    RAW_EEG = pop_loadbv(rawdata_dir, curr_setname);
    
    % add channel map
    RAW_EEG = pop_chanedit(RAW_EEG, 'lookup',chanmap_path);
    fprintf('Loading %s\n', curr_setname);

    %% re-referencing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % rereference to the average reference
    REREF_EEG = pop_reref(RAW_EEG, ref, 'keepref', 'on');
    
    %% downsampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    REREF_EEG = pop_resample(REREF_EEG, srate);
    
    %% apply filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    EEG_temp = pop_eegfiltnew(REREF_EEG, 'locutoff', locutoff, 'plotfreqz',0); 
    FILTERED_EEG = pop_eegfiltnew(EEG_temp, 'hicutoff', hicutoff ,'plotfreqz',0); 
    
    EEG_temp = pop_eegfiltnew(REREF_EEG, 'locutoff', locutoff_ica, 'plotfreqz',0); 
    FILTERED_ICA_EEG = pop_eegfiltnew(EEG_temp, 'hicutoff', hicutoff_ica ,'plotfreqz',0); 
    
    %% saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pop_saveset(FILTERED_EEG, 'filename', fullfile(output_dir, output_name));
    pop_saveset(FILTERED_ICA_EEG, 'filename', fullfile(output_dir, output_name_ica));
    
    %% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isPlot
        
        pop_eegplot(RAW_EEG, 1, 1, 1); %EEG
        pop_eegplot(FILTERED_EEG, 1, 1, 1); %EEG

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