%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------
% Goal of the script :
% Preprocessing steps: ICA
% 
% Perform ICA and removed the bad components.
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

useEOG = 0; % use EOG or not

isVisualInspect = 1; % whether or not to do a visual inspection

bad_icclasses = {'Eye'}; % bad channel classes (Eye, Muscle, Heart, etc...)
thresholdMin = 0.8; % if a component is more than XX% likely to be XX component,
                    % then it's a bad component
 
% channels
scalp_channels = 1:32;

isPlot = 1; % plot or not after processing for each subject

%% -------------------- SUBJ LOOP -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSubj = 1:nSubj
    
    % input files
    curr_subjno = subjNums(iSubj); 
    curr_inputname = sprintf('S%02d_nobadchan.set', curr_subjno);
    curr_inputname_ica = sprintf('S%02d_nobadchan_ica.set', curr_subjno);
    
    % output files
    output_dir = fullfile('results', sprintf('S%02d', curr_subjno));
    curr_outputname = sprintf('S%02d_icareject.set', curr_subjno);

    %% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    INTERP_EEG = pop_loadset('filename', curr_inputname, 'filepath', output_dir);
    INTERP_ICA_EEG = pop_loadset('filename', curr_inputname_ica, 'filepath', output_dir);

    %% perform ICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ICA_EEG = pop_runica(INTERP_ICA_EEG, 'icatype', 'runica', 'extended',1, 'chanind', scalp_channels);

    ICA_EEG.icaact = (ICA_EEG.icaweights*ICA_EEG.icasphere)*ICA_EEG.data(ICA_EEG.icachansind,:);
    
    %% Identify components with high correlation with EOG channels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % recompute ICA timecourses
    clear ICA_EEG.icaact;
    ICA_EEG = eeg_checkset(ICA_EEG, 'ica');
    n_ics = size(ICA_EEG.icaact, 1);
    
    if useEOG
        % compute correlation
        n_eog = length(EOG_channels_mat);
        
        badics_eog = zeros(n_ics, 1);
        corr_ic_eog = nan(n_ics, n_eog);

        for ieogchan = 1:n_eog
            eeg = ICA_EEG.data(EOG_channels_mat(ieogchan),:);

            for icomp = 1:n_ics
                ic = ICA_EEG.icaact(icomp, :);
                corr_tmp = corrcoef(ic, eeg);
                corr_ic_eog(icomp, ieogchan) = corr_tmp(1,2);
            end
        end

        badics_eog = any(abs(corr_ic_eog) > 0.7, 2);
        fprintf('Found %d bad ICs based on EOG.\n', sum(badics_eog));
    else
        badics_eog = 0;
    end
    
    %% visually inspect the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isVisualInspect
        [ALLEEG ICA_EEG CURRENTSET] = eeg_store(ALLEEG, ICA_EEG, CURRENTSET);

        pop_eegplot(INTERP_EEG, 1, 1, 1); %EEG
        pop_eegplot(ICA_EEG, 0, 1, 1); %Components

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

    ICA_EEG.reject.gcompreject = evalin('base', 'ALLEEG(end).reject.gcompreject');

    badics_num = strsplit(input('Bad components: ','s'));
    
    %% Identify components based on labels %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get labels
    ICA_EEG = pop_iclabel(ICA_EEG, 'default');

    % set threshold
    nClasses = length(ICA_EEG.etc.ic_classification.ICLabel.classes);
    threshold = zeros(nClasses, 2);

    bad_classes = ismember(ICA_EEG.etc.ic_classification.ICLabel.classes, bad_icclasses);

    threshold(bad_classes, 1) = thresholdMin;
    threshold(bad_classes, 2) = 1;

    % flag component
    ICA_EEG = pop_icflag(ICA_EEG, threshold);
    badics_label = ICA_EEG.reject.gcompreject;
    
    % flag components from visual inspection
    for i = 1:size(badics_num,2)
            badics_label(str2num(badics_num{i})) = 1;
    end
        
    fprintf('Found %i bad components belonging with p >= %.2f \nto classes: %s.\n',...
        sum(badics_label), ...
        thresholdMin, ...
        strjoin(ICA_EEG.etc.ic_classification.ICLabel.classes(bad_classes), ', '));

    %% Flag bad components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % update flagged ics
    flag_ics = [badics_label];
    ICA_EEG.reject.gcompreject = flag_ics;
    INTERP_EEG.reject.gcompreject = flag_ics;

    %% Remove bad components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    remove_ics = find(INTERP_EEG.reject.gcompreject);

    fprintf('Removing %d components:\n', length(remove_ics));
    fprintf(' %g', remove_ics);
    fprintf('.\n');

    ICA_EEG = pop_subcomp(ICA_EEG, remove_ics, 0);
    
    %% saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pop_saveset(ICA_EEG, 'filename', fullfile(output_dir, curr_outputname));
    
    %% documentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ica_filename = sprintf('ica_%02d.txt', curr_subjno);
    ica_file = fopen(fullfile(output_dir, ica_filename), 'w');

    fprintf(ica_file, ...
        ['For subject %2d, our ICA decomposition yields %i components. ', ...
        'From those, we rejected a total of %i components, ',...
        'with %i being eyeblinks as identified by the eeglab-generated labels, ',...
        'and %i being eye movements as identified with EOG.'],...
        curr_subjno, n_ics, length(remove_ics), sum(badics_label), sum(badics_eog));

    fclose('all');
    
    %% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isPlot
        
        pop_eegplot(INTERP_EEG, 1, 1, 1); %EEG
        pop_eegplot(ICA_EEG, 1, 1, 1); %EEG

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
