%% Brief tutorial on how to use the get_entropy() EEGLAB plugin

% launch eeglab
eeglab

% if you haven'installed the plugin yet, either go to File > Manage EEGLAB
% extensions > search get_entropy > Install
% 
% or clone the github directory in the EEGLAB plugins folder
% or donwload the github repo and unzip it in the EEGLAB plugins folder

% Load sample EEG data (4-channel wearable EEG data)
pluginPath = fileparts(which('eegplugin_entropy.m'));
EEG = pop_loadset('filename','sample_data1.set','filepath',fullfile(pluginPath,'tutorial'));

% downsample to 128 Hz to save time for tuto
EEG = pop_resample(EEG, 128);

% Use GUI to selec all parameters manually
entropy = get_entropy(EEG);

% Compute Fuzzy entropy with command line using default parameters
entropy = get_entropy(EEG,'Fuzzy entropy');

% same but only on TP channels and using variance instead of default sd
entropy = get_entropy(EEG,'Fuzzy entropy',{'TP9' 'TP10'},[],[],'Variance');

% Multiscale fuzzy entropy with default parameters, and add the 'scales' output to
% see the frequency bounds of each scale factor
[entropy, scales] = get_entropy(EEG,'Multiscale fuzzy entropy');

% same but only 15 time scales and turning off plotting (if
[entropy, scales] = get_entropy(EEG,'Multiscale fuzzy entropy',[],[],[],[],15,[],[],1);

% [~, ~, ~,MSE_sd_filt,~,scales] = get_entropy(EEG, 'Multiscale entropy', {EEG.chanlocs.labels}, 1, 2, 'Standard deviation',15,1,[],1);
