%% Run on EEG sample dataset 

clear; close all; clc
data_path = 'C:\Users\CedricCannard\Documents\eeg_data\biosemi';
cd(data_path)
eeglab; close;
locPath = fileparts(which('dipfitdefs.m'));  % for electrodes locations

filenames = dir('*.bdf');
filenames = {filenames.name}';
nFiles = length(filenames);

ApEn = nan(64,nFiles);
SampEn = nan(64,nFiles);
FuzEn = nan(64,nFiles);
FracDim = nan(64,nFiles);
PWR = nan(64, 197, nFiles);
for iFile = 1:nFiles
    % Load BDF data
    EEG = pop_biosig(fullfile(data_path, filenames{iFile}));

    % Remove auxiliary channels
    EEG = pop_select(EEG, 'rmchannel',{'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8','GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp'});

    % Import channel locations
    EEG = pop_chanedit(EEG,'lookup',fullfile(locPath,'standard_BEM','elec','standard_1005.elc'));

    % Notch filter to remove power line noise
    EEG = pop_eegfiltnew(EEG, 'locutoff',58,'hicutoff',62,'revfilt',1);

    % Highpass filter
    EEG = pop_eegfiltnew(EEG, 'locutoff',1);

    % Re-reference EEG data to infinity
    EEG = ref_infinity(EEG);

    % PSD whole range
    % [pwr, pwr_norm, psd, psd_norm, f] = compute_pwr(EEG.data, EEG.srate, .5, [1 EEG.srate/2], 2, 0);
    % figure('color', 'w','Toolbar','none','Menu','none','NumberTitle','Off');   
    % plot(f, psd_norm, 'LineWidth', 1);
    % box on; axis tight; ylabel('Power (dB)'); xlabel("Frequency (Hz)")
    % title('PSD raw'); legend({EEG.chanlocs.labels}); grid on
    % set(findall(gcf,'type','axes'),'fontSize',9,'fontweight','bold');

    % Downsample
    EEG = pop_resample(EEG,125);

    % Highpass filter
    EEG = pop_eegfiltnew(EEG, 'hicutoff',50);

    % Check power spectral density of raw data
    [pwr, pwr_norm, psd, psd_norm, f] = compute_pwr(EEG.data, EEG.srate, .5, [1 EEG.srate/2], 2, 0);
    figure('color', 'w','Toolbar','none','Menu','none','NumberTitle','Off');   
    plot(f, psd_norm, 'LineWidth', 1);
    box on; axis tight; ylabel('Power (dB)'); xlabel("Frequency (Hz)")
    title('PSD raw'); legend({EEG.chanlocs.labels}); grid on
    set(findall(gcf,'type','axes'),'fontSize',9,'fontweight','bold');
    PWR(:,:,iFile) = psd_norm;

    % Simple Heatmap of power spectrum
    figure('color', 'w');
    imagesc(f, 1:length(f), psd_norm);
    colormap("parula")
    c = colorbar; ylabel(c, 'Power (dB)','FontWeight','bold','FontSize',11,'Rotation',-90)
    xlabel('Frequency (Hz)'); ylabel('EEG channels'); title('Power spectrum heatmap')
    set(findall(gcf,'type','axes'),'fontSize',9,'fontweight','bold');

    % Correlatin matrix
    correlationMatrix = corr(psd_norm, 'Type', 'Pearson');
    figure('color', 'w');   
    imagesc(f, f, abs(correlationMatrix)); 
    colormap("parula"); 
    c = colorbar; ylabel(c, 'R^2','FontWeight','bold','FontSize',9)
    xlim([0 30]); ylim([0 30]); axis tight
    xlabel('Frequency (Hz)'); ylabel('Frequency (Hz)'); title('Correlation matrix (frequency domain)')
    set(findall(gcf,'type','axes'),'fontSize',9,'fontweight','bold');

    % Correlatin matrix
    correlationMatrix = corr(EEG.data, 'Type', 'Pearson');
    figure('color', 'w');   
    imagesc(f, f, abs(correlationMatrix)); 
    colormap("parula"); c = colorbar; ylabel(c, 'R^2','FontWeight','bold','FontSize',9)
    % xlim([0 30]); ylim([0 30]); 
    axis tight
    xlabel('EEG channels'); ylabel('EEG channels'); 
    title('Correlation matrix (time domain)')
    set(findall(gcf,'type','axes'),'fontSize',9,'fontweight','bold');

    % % Waterfall spectrogram
    % figure('color', 'w');   
    % win_length = 2 * EEG.srate; % 2-second window (500 samples)
    % noverlap = win_length / 2; % 50% overlap (250 samples)
    % nfft = win_length; % FFT size
    % for iChan = 1:EEG.nbchan
    %     nexttile
    %     spectrogram(EEG.data(iChan,:), kaiser(win_length, 5), noverlap, nfft, EEG.srate);
    %     view(-30,50); colormap parula; % parula, hot, bone
    %     xlim([1 30]); ylim([0 EEG.xmax/60])
    %     title(EEG.chanlocs(iChan).labels)
    % end
    % set(findall(gcf,'type','axes'),'fontSize',9,'fontweight','bold');


    % ApEn
    EEG = get_entropy(EEG,'Approximate entropy',[],[],[],[],[],[],[],false);
    ApEn(:,iFile) = EEG.entropy;

    % SampEn
    EEG = get_entropy(EEG,'Sample entropy',[],[],[],[],[],[],[],false);
    SampEn(:,iFile) = EEG.entropy;

    % FuzEn
    EEG = get_entropy(EEG,'Fuzzy entropy',[],[],[],[],[],[],[],false);
    FuzEn(:,iFile) = EEG.entropy;

    % Fractal votality
    EEG = get_entropy(EEG,'Fractal votality',[],[],[],[],[],[],[],false);
    FracDim(:,iFile) = EEG.entropy;
    % plot_entropy(EEG.entropy, EEG.chanlocs, 'Fractal votality', []);

    % Multiscale
    EEG = get_entropy(EEG,'Multiscale entropy',[],[],[],[],20,[],[],false);
    mulitscaleEn(:,:,iFile) = EEG.entropy;
    plot_entropy(EEG.entropy, EEG.chanlocs, 'Multiscale Entropy', []);

end




