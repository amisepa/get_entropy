%% Brief tutorial on how to use the get_entropy() EEGLAB plugin

clear; close all; clc

% launch eeglab
eeglab; close;
pop_editoptions('option_parallel', 1); % turn parrallel computing on (1) or off (0)

% if you haven'installed the plugin yet, either go to File > Manage EEGLAB
% extensions > search get_entropy > Install
% 
% or clone the github directory in the EEGLAB plugins folder
% or donwload the github repo and unzip it in the EEGLAB plugins folder

% Load provided sample EEG data from the tutorial directory 
% (several minutes of mind wandering, 64-channel Biosemi):
pluginPath = fileparts(which('eegplugin_entropy.m'));
EEG = pop_loadset('filename','sample_data_clean.set','filepath',fullfile(pluginPath,'tutorial'));
% EEG = pop_resample(EEG, 128); % downsample to 128 Hz to increase speed
% EEG = pop_select(EEG, 'point', [1 23041]);
EEG = ref_infinity(EEG);

% Launch GUI to selec all parameters manually
EEG = get_entropy(EEG);  % or Tools > Compute entropy

% Compute Fuzzy entropy with command line using default parameters
disp('-------- Fuzzy entropy ------------')
t = tic;
EEG = get_entropy(EEG,'Fuzzy entropy');
toc(t)

% If you forgot to plot outputs and want to see, you can call the function
% like this: 
plot_entropy(EEG.entropy, EEG.chanlocs, 1:EEG.nbchan)

% same but only on TP channels and using variance instead of default sd
% EEG = get_entropy(EEG,'Fuzzy entropy',{'TP9' 'TP10'},[],[],'Variance');

% Multiscale fuzzy entropy with default parameters, and add the 'scales' output to
% see the frequency bounds of each scale factor
 
% same but only 10 time scales and turning off plotting
% EEG = get_entropy(EEG,'Multiscale fuzzy entropy',[],[],[],[],10,[],[],false);
EEG = get_entropy(EEG,'Multiscale fuzzy entropy',[],[],[],[],50,[],[],false);

% same but only on channels F3 and F4 (and plotting On)
EEG = get_entropy(EEG,'Multiscale fuzzy entropy',{'F3' 'F4'},[],[],[],10,[],[],true,false);

% Same but using bandpass filtering at each time scale to control for
% spectral contamination (and plotting On)
EEG = get_entropy(EEG,'Multiscale fuzzy entropy',{'F3' 'F4'},[],[],[],10,true,[],true);

% EEG = get_entropy(EEG);   % GUI mode
% EEG = get_entropy(EEG,'Approximate entropy');
% EEG = get_entropy(EEG,'Sample entropy');
% EEG = get_entropy(EEG,'Fuzzy entropy',{'Cz'});
% EEG = get_entropy(EEG,'Multiscale entropy', {'Cz' 'O1' 'Fz'}, [],[],[],30);
% EEG = get_entropy(EEG,'Multiscale fuzzy entropy', {'Cz' 'O1' 'Fz'}, [],[],[],30);
% EEG = get_entropy(EEG, 'Sample entropy',[], 1, 2, 'Mean', 1, 2);
% EEG = get_entropy(EEG, 'Multiscale entropy', {EEG.chanlocs.labels}, 1, 2, 'Standard deviation',15,1,[],1);

