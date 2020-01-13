%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2019-12-29
% All rights reserved.
%========================================


%=========================================================================
% initialize
clear
clc
figHandles = findall(0,'Type','figure');
for n = 1:length(figHandles)
    clf(figHandles(n)) % clear figures, but don't close them
end


%=========================================================================
% settings  <== adjust based on the analysis being done
metadataFilename = 'C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2018-02-21_invivo\2018-02-21_invivo_processing.xlsx';
% metadataFilename = 'C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2018-01-31_in-vivo\2018-01-31_invivo_processing.xlsx';
metadataFilename = 'C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2018-02-21_invivo\2018-02-21_invivo_processing.xlsx';

PROCESS_RAW_DATA = 0;  % process raw data from source files
PLOT_RAW = 1;  % turn on plots when processing raw data
CREATE_SWEEPS = 1;  % create sweeps
USE_SEM = 0;  % use standard error of mean; else use standard deviation
CREATE_TIME_BACKSCATTER = 0;
LINENUM_SINGLE = 17;  % for good in-vivo backscatter (not through-tissue)

analysis_id = 'b024';
dataformat = 40;  % 10 is summer '17, 20 is jan'18, 30 is feb '18, 40 is jan'19

emgWindowProcess = [0, 40e-3];  % s
emgWindowArtifact = [0, 1.2e-3];  % s
emgWindowAmplitude = [1.2e-3, 7.4e-3];  % s
emgWindowView = [0, 9.5e-3];  % s
%emgSamplingFreq = 100000;  % Hz


mkdir(['emg_process_dataout\', analysis_id])
figuresFoldername = ['emg_process_dataout\', analysis_id, '\'];


addpath('C:\Users\david\dkpStore\UCBSF\R\Writing\Resources\figure_colorschemes\DrosteEffect-CubeHelix-00335f8\DrosteEffect-CubeHelix-00335f8')
myColormap = cubehelix(1024, 1.72, 1.31, 1.95, 1.30, [0.15, 0.82], [0.15, 0.82]);

numColors = 5;
for n = 1:numColors
    cmapIndex = (floor((n - 1).*((length(myColormap)-1)./(numColors-1)))+1);
    color5{n} = myColormap(cmapIndex, :);
%     figure(31)
%     plot((1+n):(10+n), 'Color', color5{n}, 'LineWidth', 1.5); hold on
end

set(0,'defaultAxesFontName', 'Arial');
set(0,'defaultTextFontName', 'Arial');
set(groot,'defaultAxesFontName', 'Arial');
set(groot,'defaultTextFontName', 'Arial');

%=========================================================================
% get metadata
[metadata, metadata_text, metadata_raw] = xlsread(metadataFilename);
[~,~,settingsData] = xlsread(metadataFilename, 2);
scopeFolderAndPrefix = settingsData{2,2};
emgFolderAndPrefix = settingsData{3,2};
backscatterFolderAndPrefix = settingsData{6,2};
emgSamplingFreq = settingsData{4,2};


%=========================================================================
% prep
% raw data band-pass filter
firNumTaps = 456;
lowerCutoffFreq = 200;  % Hz
upperCutoffFreq = 10000;   % Hz
nyquistRate = emgSamplingFreq./2;  
Wn = [lowerCutoffFreq./nyquistRate];
firCoefsRawHighPass = fir1(firNumTaps, Wn, 'high');
% freqz(firCoefsRawHighPass, 1, 1000);
% pause

% raw data low-pass filter
firNumTaps = 690;
cutoffFreq = 28000;   % Hz
nyquistRate = emgSamplingFreq./2;  
Wn = [cutoffFreq./nyquistRate];
firCoefsRawLowPass = fir1(firNumTaps, Wn, 'low');
% freqz(firCoefsRawLowPass, 1, 1000);
% pause

decimationFactor = 1;
emgSamplingFreqDecimated = emgSamplingFreq./decimationFactor;

% raw data high-pass filter
firNumTaps = 198;
cutoffFreq = 50;   % Hz
nyquistRate = emgSamplingFreqDecimated./2;  
Wn = [cutoffFreq./nyquistRate];
firCoefsRawHighPass = fir1(firNumTaps, Wn, 'high');
% freqz(firCoefsRawHighPass, 1, 1000);
% pause



scopeSamplingFreqOld1 = 0;
scopeSamplingFreqOld2 = 0;



%****************************************************************************************************************************
%****************************************************************************************************************************
%=========================================================================
% process each lineNum
if(PROCESS_RAW_DATA)
    data = struct('lineNum', {});
    for lineNum = 9:10  %17:28  % 8:65       %  lineNum = 4:size(metadata,1);  % [LINENUM_SINGLE]   % <== vary this for the particular analysis
        fprintf('\nlineNum: %u  ', lineNum)
        if(lineNum ~= metadata(lineNum, 1))
            fprintf('warning: lineNum %u does not match with spreadsheet; pausing\n', lineNum)
            pause
        end

        %======== get metadata
        data(lineNum).lineNum = metadata(lineNum, 1);
        data(lineNum).backscatterNum = metadata(lineNum, 2);
        data(lineNum).emgNum = metadata(lineNum, 3);
        data(lineNum).scopeNum = metadata(lineNum, 4);
        data(lineNum).tdc = metadata(lineNum, 5) .* 1e-6;  % data source in us; here in s
        data(lineNum).usDuration = metadata(lineNum, 6) .* 1e-6;  % data source in us; here in s; this is ultrasound protocol stim pulse duration; actual stim pulse is about 9 us shorter
        data(lineNum).emgGain = metadata(lineNum, 12);
        data(lineNum).vStimMax = metadata(lineNum, 9);
        data(lineNum).current = metadata(lineNum, 11) .* 1e-6;  % data source in uA; here in A
        data(lineNum).emgThreshFiltered = metadata(lineNum, 18);
        data(lineNum).emgThreshUnfiltered = metadata(lineNum, 19);
        data(lineNum).prfEst = metadata(lineNum, 8);
        
        false_triggers = [];
        eval(metadata_text{lineNum,20})  % assigns sweepLines
        data(lineNum).false_triggers = false_triggers;
        
        
        

        data(lineNum).duration = data(lineNum).usDuration - 8e-6;  % actual stim pulse is about 9 us shorter than ultrasound stim duration
        data(lineNum).chargePerPhase = data(lineNum).current .* data(lineNum).duration;  

        %======== process emg 
        if(~isnan(metadata(lineNum, 3))) % check to make sure this trial has an emg capture     

            %======== load emg data
            %data(lineNum).dataEmgRaw = csvread([emgFolderAndPrefix num2str(data(lineNum).emgNum) '.txt']);
            data(lineNum).dataEmgRaw = csvread([emgFolderAndPrefix sprintf('%02d', data(lineNum).emgNum) '.txt']);
            data(lineNum).dataEmgRaw(:,2) = data(lineNum).dataEmgRaw(:,2)./data(lineNum).emgGain;
            data(lineNum).dataEmgRaw(:,3) = data(lineNum).dataEmgRaw(:,2)./data(lineNum).emgGain;
            % data(lineNum).emgTimeStep = data(lineNum).dataEmgRaw(2,1) - data(lineNum).dataEmgRaw(1,1);
            data(lineNum).emgTimeStep = 1./emgSamplingFreq;  % s

            if(PLOT_RAW)
                figure(1); clf; hold on;
                plot(data(lineNum).dataEmgRaw(:,1), data(lineNum).dataEmgRaw(:,2));
                hold on
                plot(data(lineNum).dataEmgRaw(:,1), data(lineNum).dataEmgRaw(:,3));
                title(sprintf('EMG raw recording line %u', lineNum), 'Interpreter', 'none')
                ax = get(gca); ax.XAxis.Exponent = -3; ax.YAxis.Exponent = -3;
                xlabel('time(s)'); ylabel('EMG (V)');
                yLimCurr = ylim;
                ylim([-max(abs(yLimCurr)), max(abs(yLimCurr))])
            end
            
            %======== split into individual emg stim events and filter each
            emg_recording_style = 'continuous_no_trigger';  % 'continuous_no_trigger', 'discontinuous_trigger', 'continuous_trigger'
            if strcmp(emg_recording_style, 'continuous_no_trigger')

                firNumTaps = 698;
                cutoffFreq = 250;   % Hz
                nyquistRate = emgSamplingFreqDecimated./2;  
                Wn = [cutoffFreq./nyquistRate];
                firCoefsRawHighPassFindTrigger = fir1(firNumTaps, Wn, 'high');

                firNumTaps = 80;
                cutoffFreq = 2550;   % Hz
                nyquistRate = emgSamplingFreqDecimated./2;  
                Wn = [cutoffFreq./nyquistRate];
                firCoefsRawLowPassFindTrigger = fir1(firNumTaps, Wn, 'low');
                
                ch_to_use = 1;
                dataEmgFiltered = filtfilt(firCoefsRawHighPassFindTrigger, 1, data(lineNum).dataEmgRaw(:, (ch_to_use + 1)));
                %dataEmg = filtfilt(firCoefsRawLowPass, 1, dataEmg);
                if(lineNum == 29)
                    dataEmgFiltered(floor(82.23 .* emgSamplingFreq):floor(82.25 .* emgSamplingFreq)) = 0;   % cut out spike of noise which throws off detection of trigger points
                elseif(lineNum == 31)
                    dataEmgFiltered(floor(61.106 .* emgSamplingFreq):floor(61.112 .* emgSamplingFreq)) = 0;   % cut out spike of noise which throws off detection of trigger points
                elseif(lineNum == 32)
                    dataEmgFiltered(floor(56.88 .* emgSamplingFreq):floor(56.91 .* emgSamplingFreq)) = 0;   % cut out spike of noise which throws off detection of trigger points
                end
                dataEmgFiltered = dataEmgFiltered - mean(dataEmgFiltered);

                figure
                emgtime = data(lineNum).emgTimeStep.*(1:length(dataEmgFiltered));
                xd = data(lineNum).emgTimeStep.*(1:length(dataEmgFiltered));
                yd = data(lineNum).dataEmgRaw(:, (ch_to_use + 1));
                plot(xd(1:2:end), yd(1:2:end)); hold on
%                 plot(data(lineNum).emgTimeStep.*(1:length(dataEmgFiltered)), dataEmgFiltered)
                ax = get(gca); ax.YAxis.Exponent = -3;
                
                emg_trigger_filtered_thresh = data(lineNum).emgThreshFiltered;
                trigger_idcs_filtered = find(abs(dataEmgFiltered) > emg_trigger_filtered_thresh);
                trigger_idcs_filtered_find_diff = [0; trigger_idcs_filtered];
                trigger_idcs_diff = trigger_idcs_filtered_find_diff(2:end) - trigger_idcs_filtered_find_diff(1:(end-1));
                trigger_idcs_filtered_one_per_stim_idcs = find(trigger_idcs_diff > emgSamplingFreq * (0.1 ./ data(lineNum).prfEst));
                trigger_idcs_filtered_one_per_stim = [trigger_idcs_filtered(1); trigger_idcs_filtered(trigger_idcs_filtered_one_per_stim_idcs)];
                
%                 plot(emgtime(trigger_idcs_filtered_one_per_stim), zeros(length(trigger_idcs_filtered_one_per_stim)), 'g*')
                
                % isolate periodic triggers
                time_points_trig = emgtime(trigger_idcs_filtered_one_per_stim);
                phases_to_check = linspace(0, 2*pi, 200);
                for phase_idx = 1:length(phases_to_check)
                    likelihood = 5e-3 .* sin(2*pi*data(lineNum).prfEst*emgtime - phases_to_check(phase_idx));
                    %h1 = plot(emgtime, likelihood);
                    %pause(.001);
                    %delete(h1);
                    %mult(phase_idx) = likelihood' .* abs(dataEmgFiltered);
                    mult(phase_idx) = dot(likelihood, dataEmgFiltered.^2);
                end
                [~,phase_to_use_idx] = max(mult);
                likelihood_to_use = 1 .* sin(2*pi*data(lineNum).prfEst*emgtime - phases_to_check(phase_to_use_idx));
                
                thresh_likelihood = 0.9;
                trigger_idcs_filtered_one_per_stim(likelihood_to_use(trigger_idcs_filtered_one_per_stim) < thresh_likelihood) = [];
                
%                 plot(emgtime(trigger_idcs_filtered_one_per_stim), zeros(length(trigger_idcs_filtered_one_per_stim)), 'r*')
                
                emg_trigger_unfiltered_thresh = data(lineNum).emgThreshUnfiltered;
                for stim_pulse_idx = 1:length(trigger_idcs_filtered_one_per_stim)
                    range_search = floor(trigger_idcs_filtered_one_per_stim(stim_pulse_idx) + emgSamplingFreq .* (-0.02)):floor(trigger_idcs_filtered_one_per_stim(stim_pulse_idx) + emgSamplingFreq .* 0.02);
                    %diffs = data(lineNum).dataEmgRaw(range_search+1, (ch_to_use + 1)) - data(lineNum).dataEmgRaw(range_search+0, (ch_to_use + 1));
                    over_thresh = find(abs(   data(lineNum).dataEmgRaw(range_search, (ch_to_use + 1)) - data(lineNum).dataEmgRaw(range_search(1), (ch_to_use + 1))) > emg_trigger_unfiltered_thresh);
                    trigger_idcs_real(stim_pulse_idx) = range_search(1) + over_thresh(1);
                end
                
                plot(emgtime(trigger_idcs_real), zeros(length(trigger_idcs_real)), 'c*')
                title(num2str(lineNum))
                
                emgPulseStartIndices = trigger_idcs_real;
                
                clearvars trigger_idcs_filtered_one_per_stim    trigger_idcs_real    mult    trigger_idcs_filtered_one_per_stim_idcs;
                
                

            elseif strcmp(emg_recording_style, 'dicontinuous_trigger')
                timeDifferences = data(lineNum).dataEmgRaw(2:end,1) - data(lineNum).dataEmgRaw(1:(end-1),1);
        %         figure(2)
        %         plot(timeDifferences)

                emgPulseStartIndices = find(timeDifferences > 20.*data(lineNum).emgTimeStep) + 1;
                emgPulseStartIndices = emgPulseStartIndices(2:(end));

                if(PLOT_RAW)
                    figure(1);
                    plot(data(lineNum).dataEmgRaw(emgPulseStartIndices, 1), 0, 'g*')
                end
                
            elseif strcmp(emg_recording_style, 'continuous_trigger')
                ;
            end
            
            emgWindowProcessIndices = floor(emgWindowProcess(1)./data(lineNum).emgTimeStep):1:ceil(emgWindowProcess(2)./data(lineNum).emgTimeStep);  % s
            emgWindowArtifactIndices = 1 + floor(emgWindowArtifact(1)./data(lineNum).emgTimeStep):1:floor(emgWindowArtifact(2)./data(lineNum).emgTimeStep); 
            emgWindowAmplitudeIndices = 1 + floor(emgWindowAmplitude(1)./data(lineNum).emgTimeStep):1:floor(emgWindowAmplitude(2)./data(lineNum).emgTimeStep);  % s
            emgWindowViewIndices = 1 + floor(emgWindowView(1)./data(lineNum).emgTimeStep):1:floor(emgWindowView(2)./data(lineNum).emgTimeStep);  % s

            data(lineNum).emgTimeTime = data(lineNum).emgTimeStep .* emgWindowProcessIndices;
            data(lineNum).emgTimeTimeDecimated = data(lineNum).emgTimeTime(1:decimationFactor:end);
            data(lineNum).numEmgPulses = length(emgPulseStartIndices);


            emgStartTimes = data(lineNum).dataEmgRaw(emgPulseStartIndices, 1);
            emgStartTimeDifferences = emgStartTimes(2:end) - emgStartTimes(1:(end-1));
            data(lineNum).prf = 1./median(emgStartTimeDifferences);

            if(PLOT_RAW); figure(3); clf; end
            
            dataEmg_LPF = filtfilt(firCoefsRawLowPassFindTrigger, 1, data(lineNum).dataEmgRaw(:, (ch_to_use + 1)));

            
            for n = 1:data(lineNum).numEmgPulses
    %             dataEmg = (data(lineNum).dataEmgRaw(emgPulseStartIndices(n):(emgPulseStartIndices(n)+numelEmgPulse-1), 2))';
                dataEmg = (data(lineNum).dataEmgRaw((emgPulseStartIndices(n) + emgWindowProcessIndices), 2))';
                %dataEmg = dataEmg - mean(dataEmg);
                %dataEmg = dataEmg - dataEmg(1);
                dataEmg = dataEmg - mean(data(lineNum).dataEmgRaw((emgPulseStartIndices(n) + (floor(emgSamplingFreq .* -0.005):floor(emgSamplingFreq .* -0.0005))), 2));
    %             subplot(2,1,1); plot(data(lineNum).emgTimeTime, dataEmg); hold on;
    %             dataEmg = filtfilt(firCoefsRawHighPass, 1, dataEmg);
    %             dataEmg = filtfilt(firCoefsRawLowPass, 1, dataEmg);
                dataEmg = dataEmg(1:decimationFactor:end);
    %             dataEmg = filtfilt(firCoefsRawHighPass, 1, dataEmg);
    %             dataEmg = dataEmg - mean(dataEmg);
                dataEmg = dataEmg - dataEmg(1);
                data(lineNum).dataEmg(n,:) = dataEmg;
    %             subplot(2,1,2); 
                if(PLOT_RAW)
                    plot(data(lineNum).emgTimeTimeDecimated, dataEmg); hold on;
                end
            end
            yLimCurr = ylim;
            ylim([-max(abs(yLimCurr)), max(abs(yLimCurr))])
            ax = get(gca); ax.XAxis.Exponent = -3; ax.YAxis.Exponent = -3;  
            
            fn_format_and_save_figure(3, [pwd filesep figuresFoldername filesep sprintf('emgpop_line_%u', lineNum)], [0])
            


            %======== get single time-domain emg line for each condition
            data(lineNum).emgTimeMean = mean(data(lineNum).dataEmg, 1);
            data(lineNum).emgTimeStd = std(data(lineNum).dataEmg, 0, 1);
            data(lineNum).emgTimeSem = data(lineNum).emgTimeStd ./ sqrt(size(data(lineNum).dataEmg, 1));

            if(PLOT_RAW)
                figure(4); hold off;
                shadedErrorBar(data(lineNum).emgTimeTimeDecimated, data(lineNum).emgTimeMean, data(lineNum).emgTimeSem);
                ax = get(gca); ax.XAxis.Exponent = -3; ax.YAxis.Exponent = -3;
                title(sprintf('EMG population recording'), 'Interpreter', 'none')
                xlim(emgWindowView)
                xlabel('time(s)'); ylabel('EMG (V)');
                yLimCurr = ylim;
                ylim([-max(abs(yLimCurr)), max(abs(yLimCurr))])
            end

            %======== get artifact and emg amplitude stats of each emg
            artifactROI = data(lineNum).dataEmg(:, emgWindowArtifactIndices);
            emgROI = data(lineNum).dataEmg(:, emgWindowAmplitudeIndices);

            % artifact peak-peak
            data(lineNum).artifactAmplitudes = max(artifactROI, [], 2) - min(artifactROI, [], 2);
            data(lineNum).artifactAmplitudeMean = mean(data(lineNum).artifactAmplitudes);
            data(lineNum).artifactAmplitudeStd = std(data(lineNum).artifactAmplitudes);
            data(lineNum).artifactAmplitudeSem = data(lineNum).artifactAmplitudeStd ./ sqrt(data(lineNum).numEmgPulses);

            % emg amplitude peak=peak
            data(lineNum).emgAmplitudes = max(emgROI, [], 2) - min(emgROI, [], 2);
            data(lineNum).emgAmplitudeMean = mean(data(lineNum).emgAmplitudes);
            data(lineNum).emgAmplitudeStd = std(data(lineNum).emgAmplitudes);
            data(lineNum).emgAmplitudeSem = data(lineNum).emgAmplitudeStd ./ sqrt(data(lineNum).numEmgPulses);
            
            % emg amplitude single-sided (just above baseline)
            data(lineNum).emgAmplitudesSingle = max(abs(emgROI), [], 2);   % ********** choose max or min based on which way the primary emg deflection is away from zero
            data(lineNum).emgAmplitudeSingleMean = mean(data(lineNum).emgAmplitudesSingle);
            data(lineNum).emgAmplitudeSingleStd = std(data(lineNum).emgAmplitudesSingle);
            data(lineNum).emgAmplitudeSingleSem = data(lineNum).emgAmplitudeSingleStd ./ sqrt(data(lineNum).numEmgPulses);
            
            % emg area peak peak
            data(lineNum).emgAreas = sum(abs(emgROI), 2).*data(lineNum).emgTimeStep;  % in V s
            data(lineNum).emgAreaMean = mean(data(lineNum).emgAreas);
            data(lineNum).emgAreaStd = std(data(lineNum).emgAreas);
            data(lineNum).emgAreaSem = data(lineNum).emgAreaStd ./ sqrt(data(lineNum).numEmgPulses);
            
            % emg area single-sided (just above baseline)
            data(lineNum).emgAreasSingle = sum(emgROI(emgROI > 0)).*data(lineNum).emgTimeStep;  % in V s
            data(lineNum).emgAreaMean = mean(data(lineNum).emgAreasSingle);
            data(lineNum).emgAreaStd = std(data(lineNum).emgAreasSingle);
            data(lineNum).emgAreaSem = data(lineNum).emgAreaStd ./ sqrt(data(lineNum).numEmgPulses);

%             figure(1); saveas(gcf, sprintf('quality_check_figures/all_%u.png', lineNum));
%             figure(3); saveas(gcf, sprintf('quality_check_figures/overlaid_%u.png', lineNum));
%             figure(4); saveas(gcf, sprintf('quality_check_figures/timeErrorBar_%u.png', lineNum));
            data(lineNum).dataEmgRaw = [];  % clear dataEmgRaw from structure to save space when saving
        end

        
        %======== process scope data
        if(~isnan(metadata(lineNum, 4)))  % check to make sure this trial has a scope capture
            scopeData = csvread([scopeFolderAndPrefix num2str(data(lineNum).scopeNum) '.csv'], 2, 0);
            
            scopeData = scopeData(1:(end-50), :);  % cut off a bit of data at the end that are zeros
            
            %vertical align
            vstimTime = scopeData(:,1);
            scopeSampleTimeStep = vstimTime(2) - vstimTime(1);
            vstim = scopeData(:,4);
            vstim = vstim - median(vstim(1:1000));
            
            vpz1 = scopeData(:,2);
            vpz2 = scopeData(:,3);
            vvdd = scopeData(:,5);
            
            % temporal align
            triggerThreshold = .175;
            zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);  
            zx = zci(vstim - triggerThreshold);
            vstimTime = vstimTime - vstimTime(zx(1));
            

            % create vstim and vdd low-pass filter and apply filter
            firNumTaps = 900;
            cutoffFreq = 2400000;   % Hz
            scopeSamplingFreq = 1./scopeSampleTimeStep;
            nyquistRate = scopeSamplingFreq./2;  
            Wn = [cutoffFreq./nyquistRate];
            if(scopeSamplingFreq ~= scopeSamplingFreqOld1)
%                 firCoefsScopeLowPass1 = fir1(firNumTaps, Wn, 'low', blackmanharris(firNumTaps+1));
%                 firCoefsScopeLowPass_new = firCoefsScopeLowPass1;
%                 negRatio = .33;
%                 firCoefsScopeLowPass_new(firCoefsScopeLowPass_new < 0) = firCoefsScopeLowPass_new(firCoefsScopeLowPass_new < 0) * negRatio;\
%                 firCoefsScopeLowPass_new = firCoefsScopeLowPass_new ./ sum(firCoefsScopeLowPass_new);
%                 firCoefsScopeLowPass1 = firCoefsScopeLowPass_new;
        
                firCoefsScopeLowPass1 = fn_gaussfiltcoef(scopeSamplingFreq, cutoffFreq);

            end            
            scopeSamplingFreqOld1 = scopeSamplingFreq;
            vstim = filtfilt(firCoefsScopeLowPass1, 1, vstim);
            vvdd = filtfilt(firCoefsScopeLowPass1, 1, vvdd);
            
            % create pz1 pz2 low-pass filter and apply filter
            firNumTaps = 900;
            cutoffFreq = 8400000;   % Hz
            scopeSamplingFreq = 1./scopeSampleTimeStep;
            nyquistRate = scopeSamplingFreq./2;  
            Wn = [cutoffFreq./nyquistRate];
            if(scopeSamplingFreq ~= scopeSamplingFreqOld2)
%                 firCoefsScopeLowPass2 = fir1(firNumTaps, Wn, 'low');
                firCoefsScopeLowPass2 = fn_gaussfiltcoef(scopeSamplingFreq, cutoffFreq);
            end     
            scopeSamplingFreqOld2 = scopeSamplingFreq;
            vpz1 = filtfilt(firCoefsScopeLowPass2, 1, vpz1);
            vpz2 = filtfilt(firCoefsScopeLowPass2, 1, vpz2);

            % plot
            if(PLOT_RAW)
            figure(6); hold off;
            plot(vstimTime, vstim, 'DisplayName', sprintf('tdc %u, lineNum %u', data(lineNum).tdc .* 1e6, data(lineNum).lineNum))
%             xlim([-.015e-3, .479e-3])
            lgd = legend('show', 'boxoff', 'Location', 'northeast');
            lgd.FontSize = 10;
            ax = get(gca); ax.XAxis.Exponent = -3; ax.YAxis.Exponent = 0;
%             figure(6); saveas(gcf, sprintf('quality_check_figures/vstim_%u.png', lineNum));
            end
            
            % input vstim data to data struct
            data(lineNum).vstimTime = vstimTime;
            data(lineNum).vstim = vstim;
            data(lineNum).vpz1 = vpz1;
            data(lineNum).vpz2 = vpz2;
            data(lineNum).vvdd = vvdd;
            
        end
        
        %======== process backscatter data
        if(~isnan(data(lineNum).backscatterNum))  % check to make sure this trial has a scope capture
            backscatterFilenames = dir([backscatterFolderAndPrefix, '*']);
            backscatterFilenames = {backscatterFilenames(:).name};
            tf = contains(backscatterFilenames, ['cond', sprintf('%02d', data(lineNum).backscatterNum)]);
            bfn_to_use = find(tf);
            backscatterFilename = backscatterFilenames{bfn_to_use(1)};
            
            %backscatterLoaded = load([backscatterFolderAndPrefix num2str(data(lineNum).backscatterNum) '.mat']);
            backscatterLoaded = load([backscatterFolderAndPrefix, backscatterFilename]);
            
            backscatterLoaded.samplingFreq = metadata(lineNum, 17);    % ********** only for 2018-01-31
            
            backscatterScale = (0.03./0.0896).*0.8/4096.0;   % 0.8 V range divided into 2^12 steps with front-end gain of .0896/.03
    
            for n = 1:length(backscatterLoaded.rawSignalStoreB)
                if n == length(backscatterLoaded.rawSignalStoreB) - 1
                    plotOn = 1;
                else
                    plotOn = 0;
                end
                
                rawBackscatterScaled = backscatterScale .* (backscatterLoaded.rawSignalStoreB(n).rawSignal);
                
                [out1, out2, out3, out4, out5] = fn_stimBackscatterProcess(rawBackscatterScaled', ...
                                                                        backscatterLoaded.samplingFreq, ...
                                                                        1.8e6, ...
                                                                        data(lineNum).tdc, ...
                                                                        data(lineNum).usDuration, ...
                                                                        [8.5e-6, 35e-6], ...   % 9, 34
                                                                        plotOn);
                if(n==1)
                    data(lineNum).backscatterTime = out1;
                    data(lineNum).backscatterSamplingFreq = backscatterLoaded.samplingFreq;
                    currRow = 1;
                end
                if(length(out1) == length(data(lineNum).backscatterTime))
                    data(lineNum).backscatter(currRow,:) = out2;
                    data(lineNum).backscatterL1normRecharge(currRow) = out3;
                    data(lineNum).backscatterL1normStim(currRow) = out4;
                    data(lineNum).backscatterL1normRatio(currRow) = out5;
                    currRow = currRow + 1;
                end
                if n == length(backscatterLoaded.rawSignalStoreB) - 1
%                     figure(21); saveas(gcf, sprintf('quality_check_figures/backscatter2_%u.png', lineNum));
    %               figure(22); saveas(gcf, sprintf('quality_check_figures/backscatterFFT_%u.png', lineNum));
                end    
            end
        end
    end
%     save('dataWorkingThroughTissue.mat', 'data');
%     save('dataWorkingBackscatter1.mat', 'data');
%     save('dataWorkingThroughTissue3.mat', 'data');
    %save('emg_process_dataout\dataWorking.mat', 'data');
    save([figuresFoldername, 'dataWorking.mat'], 'data');
%     save('single.mat', 'data');
%     save('data150.mat', 'data');

end

















%****************************************************************************************************************************
%****************************************************************************************************************************


%=========================================================================
% figure settings
figSize = [900 550];
figFontSize = 16;
figTitleSize = 4;
figAxislabelSize = 4;
figTicklabelSize = 4;



%=========================================================================
% plot sweeps 
if(CREATE_SWEEPS) 
    analysisDataLoad = 'b023';
    load(['emg_process_dataout\', analysisDataLoad, '\dataWorking.mat'])
    eval(settingsData{5,2})  % assigns sweepLines
    
    for sweepNum = sweepLines
        eval(settingsData{sweepNum, 2})  % assigns lineNumsSweep

        CURRENT_SWEEP = 0;
        DURATION_SWEEP = 0;
        if(strcmp(settingsData{sweepNum, 3}, 'current'))
            CURRENT_SWEEP = 1;
        elseif(strcmp(settingsData{sweepNum, 3}, 'duration'))
            DURATION_SWEEP = 1;
        else
            disp('neither current nor duration sweep')
        end
        maxDuration = 0;
        figure(4); clf; figure(5); clf; figure(6); clf; figure(7); clf; figure(8); clf; figure(9); clf;
        recruitAmpPopX = [];
        recruitAmpPopY = [];
        
        %======== time-domain plot and collect data
        for n = 1:length(lineNumsSweep)
            if(length(lineNumsSweep{n}) > 1)
               ;  % merge lineNums here 
            end
            lineNum = lineNumsSweep{n};
            cmapIndex = (length(myColormap) + 1) - (floor((n - 1).*((length(myColormap)-1)./(length(lineNumsSweep)-1)))+1);
            colorThisCondition = myColormap(cmapIndex, :);

            %======== plot time-domain
            if(~isnan(metadata(lineNum, 3)))  % check to make sure this trial has a scope capture
                figure(4); hold on;
                
                fprintf('current: %d     duration:    %d      N:     %d     \n', data(lineNum).current, data(lineNum).duration, size((data(lineNum).dataEmg), 1))
                if(USE_SEM)
                    plotDataErr = data(lineNum).emgTimeSem;  %  .* 1000
                else
                    plotDataErr = data(lineNum).emgTimeStd;  %  .* 1000
                end
                
                if(CURRENT_SWEEP)
                    H = shadedErrorBar(data(lineNum).emgTimeTimeDecimated .* 1000, data(lineNum).emgTimeMean .* 1000, plotDataErr .* 1000, 'lineProps', {'-k', 'Color', colorThisCondition, 'LineWidth', 1.00, 'DisplayName', sprintf('%3.0f %cA', data(lineNum).current .* 1e6, 956)});  % , 'transparent', false
                elseif(DURATION_SWEEP)
                    H = shadedErrorBar(data(lineNum).emgTimeTimeDecimated .* 1000, data(lineNum).emgTimeMean .* 1000, plotDataErr .* 1000, 'lineProps', {'-k', 'Color', colorThisCondition, 'LineWidth', 1.00, 'DisplayName', sprintf('%3.0f %cs', data(lineNum).duration .* 1e6, 956)});
                else
                    H = shadedErrorBar(data(lineNum).emgTimeTimeDecimated .* 1000, data(lineNum).emgTimeMean .* 1000, plotDataErr .* 1000, 'lineProps', {'-k', 'Color', colorThisCondition, 'LineWidth', 1.00, 'DisplayName', sprintf('%3.0f %cA', data(lineNum).current .* 1e6, 956)});
                end

                HforLegend(n) = H.mainLine;
                ax = get(gca); ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
                title(sprintf('Compound muscle action potential \n stimulation current sweep'), 'Interpreter', 'none')
                xlim(emgWindowView .* 1000)
                xlabel('time (ms)'); ylabel('EMG (mV)');
                yLimCurr = ylim;
                ylim([-max(abs(yLimCurr)), max(abs(yLimCurr))])
                ylim([-4.5, 4.5])
            end
            
            %======== plot vStim trace
            if(~isnan(metadata(lineNum, 4)))  % check to make sure this trial has a scope capture
                % plot
                figure(6); hold on;
                if(CURRENT_SWEEP)
%                     indicesToPlot = find((data(lineNum).vstimTime))
                    indicesToPlot = (data(lineNum).vstimTime > -7e-6) & (data(lineNum).vstimTime < 474e-6);
                    hvstim = plot(data(lineNum).vstimTime(indicesToPlot) .* 1e6, data(lineNum).vstim(indicesToPlot), 'DisplayName', sprintf('%3.0f %cA', data(lineNum).current .* 1e6, 956), 'Color', colorThisCondition);
                    xlim(1e6.*[-35e-6, (data(lineNum).duration + 80e-6 + 25e-6)])
                    title('stim voltage waveform vs stim current')
                elseif(DURATION_SWEEP)
                    hvstim = plot(data(lineNum).vstimTime .* 1e6, data(lineNum).vstim, 'DisplayName', sprintf('%3.0f %cs', data(lineNum).duration .* 1e6, 956), 'Color', colorThisCondition);
                    maxDuration = max(maxDuration, data(lineNum).duration);
                    xlim(1e6.*[-60e-6, (maxDuration + 80e-6 + 60e-6)])
                    title('stim voltage waveform vs stim duration')
                else
                    hvstim = plot(data(lineNum).vstimTime .* 1e6, data(lineNum).vstim, 'DisplayName', sprintf('%3.0f %cA', data(lineNum).current .* 1e6, 956), 'Color', colorThisCondition);
                    xlim(1e6.*[-60e-6, (data(lineNum).duration + 80e-6 + 60e-6)])
                    title('stim voltage waveform vs stim current')
                end
%                 xlim([-.015e-3, .479e-3])
                ax = get(gca); ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
                set(hvstim, 'LineWidth', 1.5)
                xlabel('time (\mus)')
                ylabel('stim voltage (V)')
            end
            
            %======== collect sweep data
            vStimX(n) = data(lineNum).tdc;
            vStimY(n) = data(lineNum).vStimMax;
            
            %======== X data
            if(CURRENT_SWEEP)
                recruitArtX(n) = data(lineNum).current;
                recruitAmpX(n) = data(lineNum).current;
            elseif(DURATION_SWEEP)
                recruitArtX(n) = data(lineNum).duration;
                recruitAmpX(n) = data(lineNum).duration;
            else
                recruitArtX(n) = data(lineNum).current;
                recruitAmpX(n) = data(lineNum).current;
            end
            
            
%             recruitArtY(n) = data(lineNum).artifactAmplitudeMean .* 1000;
%             recruitAmpY(n) = data(lineNum).emgAmplitudeMean .* 1000;
%             if(USE_SEM)
%                 recruitArtErr(n) = data(lineNum).artifactAmplitudeSem .* 1000;
%                 recruitAmpErr(n) = data(lineNum).artifactAmplitudeSem .* 1000;
%             else
%                 recruitArtErr(n) = data(lineNum).artifactAmplitudeStd .* 1000;
%                 recruitAmpErr(n) = data(lineNum).artifactAmplitudeStd .* 1000;
%             end
            
            %======== Y and error data
            recruitArtY(n) = data(lineNum).artifactAmplitudeMean;
            if(USE_SEM)
                recruitArtErr(n) = data(lineNum).artifactAmplitudeSem;
            else
                recruitArtErr(n) = data(lineNum).artifactAmplitudeStd;
            end
            
            recruitAmpY(n) = data(lineNum).emgAmplitudeSingleMean;
            if(USE_SEM)
                recruitAmpErr(n) = data(lineNum).emgAmplitudeSingleSem;
            else
                recruitAmpErr(n) = data(lineNum).emgAmplitudeSingleStd;
            end
            
            
            %======== population data
            recruitAmpPopY = [recruitAmpPopY; data(lineNum).emgAmplitudesSingle];
            if(CURRENT_SWEEP)
                recruitAmpPopX = [recruitAmpPopX; (ones((length(recruitAmpPopY) - length(recruitAmpPopX)), 1).*data(lineNum).current)];
            elseif(DURATION_SWEEP)
                recruitAmpPopX = [recruitAmpPopX; (ones((length(recruitAmpPopY) - length(recruitAmpPopX)), 1).*data(lineNum).duration)];
            else
                recruitAmpPopX = [recruitAmpPopX; (ones((length(recruitAmpPopY) - length(recruitAmpPopX)), 1).*data(lineNum).duration)];
            end

        end
        
        figure(4)
%         lgd4 = legend('show', 'Location', 'northeast');
        lgd4 = legend(HforLegend);
        legend boxoff;
        lgd4.FontSize = 12; lgd4.EdgeColor = [.65, .65, .65];
        lgd4.Location = 'eastoutside';
        
        figure(6)
        lgd6 = legend('show');
        legend boxoff;
        lgd6.Location = 'eastoutside';
        lgd6.FontSize = 12; lgd6.EdgeColor = [.65, .65, .65];        

        

        
        %======== recruitment curve sigmoid fit
        [xData, yData] = prepareCurveData( recruitAmpPopX, recruitAmpPopY );

        % Set up fittype and options.
        ft = fittype( 'a/(1+exp(-b*(x-c))) + d', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.DiffMinChange = 1e-12;
        opts.Display = 'Off';
        opts.Lower = [0 -Inf -Inf 0];
        opts.MaxFunEvals = 60000;
        opts.MaxIter = 40000;
        opts.StartPoint = [0.004 50000 0.0001 1e-06];
        opts.TolFun = 1e-12;
        opts.TolX = 1e-08;

        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );

        % Plot fit with data.
%         figure( 'Name', 'untitled fit 1' );
%         h = plot( fitresult, xData, yData );
%         legend( h, 'recruitAmpPopY vs. recruitAmpPopX', 'untitled fit 1', 'Location', 'NorthEast' );
%         % Label axes
%         xlabel recruitAmpPopX
%         ylabel recruitAmpPopY
%         grid on        

        x = -1e-3:.1e-6:1e-3;
        yFit = fitresult.a./(1+exp(-fitresult.b.*(x - fitresult.c))) + fitresult.d;
        
        
        %======== plot recruitment curve
        colorLineSweeps = [0, 0.1563, 1.0000];
        
        figure(7); hold off
        errorbar(recruitArtX .* 1e6, recruitArtY .* 1e3, recruitArtErr .* 1e3, '-', 'LineWidth', 1.5, 'Color', colorLineSweeps)
        ax = get(gca); ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
        ylabel('artifact amplitude (mV)')

        figure(8); clf; hold on;
        plot(x .* 1e6, yFit .* 1e3, 'Color', color5{(end-1)}, 'LineWidth', 1.25);
        errorbar(recruitAmpX .* 1e6, recruitAmpY .* 1e3, recruitAmpErr .* 1e3, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', .1, 'CapSize', 10, 'LineWidth', 1.25)  % data in s and V; plot in us and mV
        plot(recruitAmpPopX .* 1e6, recruitAmpPopY .* 1e3, 'b.', 'MarkerEdgeColor', color5{2}, 'MarkerFaceColor', color5{2}, 'MarkerSize', 6)
        ax = get(gca); ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
        ylabel('emg amplitude (mV)')
        xlim([(min(recruitAmpX .* 1e6) - 50) (max(recruitAmpX .* 1e6) + 50)])
        ylim(1000.*[-.0006, .0045])
        
        figure(9); hold off
        H9 = shadedErrorBar(recruitAmpX .* 1e6, recruitAmpY .* 1e3, recruitAmpErr .* 1e3, 'lineProps', {'-k', 'Color', colorLineSweeps, 'DisplayName', sprintf('%3.0f %cA', data(lineNum).current .* 1e6, 956)});
        H9forLegend = H.mainLine;
        ax = get(gca); ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
        ylabel('emg amplitude (mV)')
        ylim(1000.*[-.001, .005])
        
        if(CURRENT_SWEEP)
            figure(7); xlabel('stim current (\muA)'); title('artifact amplitude vs stim current'); xlim(1e6 .* [0 450e-6])
            figure(8); xlabel('stim current (\muA)'); title('emg amplitude vs stim current'); xlim(1e6 .* [0 450e-6])
            figure(9); xlabel('stim current (\muA)'); title('emg amplitude vs stim current'); xlim(1e6 .* [0 450e-6])
        elseif(DURATION_SWEEP)
            figure(7); xlabel('stim duration (\mus)'); title('artifact amplitude vs stim duration'); xlim(1e6 .* [0 430e-6])
            figure(8); xlabel('stim duration (\mus)'); title('emg amplitude vs stim duration'); xlim(1e6 .* [0 430e-6])
            figure(9); xlabel('stim duration (\mus)'); title('emg amplitude vs stim duration'); xlim(1e6 .* [0 430e-6])
        else
            figure(7); xlabel('stim current (\muA)'); title('artifact amplitude vs stim current'); xlim(1e6 .* [0 450e-6])
            figure(8); xlabel('stim current (\muA)'); title('emg amplitude vs stim current'); xlim(1e6 .* [0 450e-6])
            figure(9); xlabel('stim current (\muA)'); title('emg amplitude vs stim current'); xlim(1e6 .* [0 450e-6])
        end
        
        
        %======== save figures
        xlabel('X Axis', 'FontSize',12)
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 16)
      
%         
%         figure(4); set(gcf, 'Color', 'White', 'pos', [40, 40, 40+figSize(1), 40+figSize(2)]); box off; drawnow; set(gca, 'FontSize', figFontSize);
%         figure(5); set(gcf, 'Color', 'White', 'pos', [40, 40, 40+figSize(1), 40+figSize(2)]); box off; drawnow; set(gca, 'FontSize', figFontSize);
%         figure(6); set(gcf, 'Color', 'White', 'pos', [40, 40, 40+figSize(1), 40+figSize(2)]); box off; drawnow; set(gca, 'FontSize', figFontSize);
%         figure(7); set(gcf, 'Color', 'White', 'pos', [40, 40, 40+figSize(1), 40+figSize(2)]); box off; drawnow; set(gca, 'FontSize', figFontSize);
%         figure(8); set(gcf, 'Color', 'White', 'pos', [40, 40, 40+figSize(1), 40+figSize(2)]); box off; drawnow; set(gca, 'FontSize', figFontSize);
%         figure(9); set(gcf, 'Color', 'White', 'pos', [40, 40, 40+figSize(1), 40+figSize(2)]); box off; drawnow; set(gca, 'FontSize', figFontSize);
        
        
%         figure(4); saveas(gcf, sprintf('figures_out1/emgpop_sweep_%u.png', sweepNum))
%         figure(6); saveas(gcf, sprintf('figures_out1/vstim_sweep_%u.png', sweepNum))
%         figure(7); saveas(gcf, [pwd filesep figuresFoldername filesep sprintf('artamp_sweep_%u.png', sweepNum)])
%         figure(8); saveas(gcf, [pwd filesep figuresFoldername filesep sprintf('emgamp_sweep_%u.png', sweepNum)])
%         figure(9); saveas(gcf, [pwd filesep figuresFoldername filesep sprintf('emgampShaded_sweep_%u.png', sweepNum)])
        
        
        fn_format_and_save_figure(4, [pwd filesep figuresFoldername filesep sprintf('emgpop_sweep_%u', sweepNum)], [0])
        fn_format_and_save_figure(6, [pwd filesep figuresFoldername filesep sprintf('vstim_sweep_%u', sweepNum)], [0])
%         fn_format_and_save_figure(7, [pwd filesep figuresFoldername filesep sprintf('artamp_sweep_%u', sweepNum)], [0])
        fn_format_and_save_figure(8, [pwd filesep figuresFoldername filesep sprintf('emgamp_sweep_%u', sweepNum)], [0])
%         fn_format_and_save_figure(9, [pwd filesep figuresFoldername filesep sprintf('emgampShaded_sweep_%u', sweepNum)], [0])
        
        clear H;
        close 7; close 9;  % close figures we don't care about

    end
       
end



















%****************************************************************************************************************************
%****************************************************************************************************************************
%=========================================================================
% plot time backscatter
 
if(CREATE_TIME_BACKSCATTER)
    
%     load('dataWorkingThroughTissue.mat');
%     load('data150.mat');
    load('dataWorking.mat');
    
    eval(settingsData{28,2})  % assigns lineNum for through-tissue
    
    
    %======= figure for overlaid time backscatter
    figure(12); hold on
    
    % note: t = 0 for this plot is beginning of stim
    % backscatter will have to be shifted back by length of time-of-flight
    % and 16-cycle hedaer
    %
    % scope values are already set to 0 at beginning of stim
    %
    % emg will have to be shifted left by length of time-of-flight and
    % 16-cycle header
    

    
    leftShiftTimeOfFlight = 16.86e-6;  % s
    leftShiftHeader = 8e-6;  % s
    
    hvstim = plot(data(lineNum).vstimTime .* 1000, data(lineNum).vstim);  % source data is ~0 to 1 V   plot in ms and V
    hvpz1 = plot(data(lineNum).vstimTime .* 1000, data(lineNum).vpz1);  % source data is ~-3 to 3 V   plot in ms and V
%     hvpz2 = plot(data(lineNum).vstimTime .* 1000, data(lineNum).vpz2);  % source data is ~-3 to 3 V   plot in ms and V
%     set(hvpz1,'Visible','off')
%     set(hvpz2,'Visible','off')    
    hvvdd = plot(data(lineNum).vstimTime .* 1000, data(lineNum).vvdd);  % source data is ~0 to 3 V   plot in ms and V
%     hvpz1vpz2 = plot(data(lineNum).vstimTime .* 1000, (data(lineNum).vpz1 - data(lineNum).vpz2));  % time data in s, plot in ms; volt data in volt, plot in volt; plot in ms and V
    
    backscatterWindow = [2e-6, 45e-6];
    freqSampling = data(lineNum).backscatterSamplingFreq;
    backscatterWindowIndicesRecharge = floor(backscatterWindow(1) .* freqSampling):ceil(backscatterWindow(2) .* freqSampling);
    backscatterWindowIndicesStim = backscatterWindowIndicesRecharge + floor((data(lineNum).tdc + data(lineNum).usDuration) .* freqSampling);
    backscatterToPlotTime = (data(lineNum).backscatterTime - (data(lineNum).tdc + leftShiftHeader + leftShiftTimeOfFlight));
    backscatterToPlotVolt = data(lineNum).backscatter(end,:);
    backscatterToPlotTime = [backscatterToPlotTime(backscatterWindowIndicesRecharge) NaN backscatterToPlotTime(backscatterWindowIndicesStim)];
    backscatterToPlotVolt = [backscatterToPlotVolt(backscatterWindowIndicesRecharge) NaN backscatterToPlotVolt(backscatterWindowIndicesStim)];
    hbackscatter = plot(backscatterToPlotTime .* 1000, backscatterToPlotVolt .* 1000 .* 0.1);  % source data is limited to ~-2048 to +2048; typically ~-800 to +800; plot in ms and mV
%     
    hemg = plot((data(lineNum).emgTimeTimeDecimated  - (leftShiftTimeOfFlight + leftShiftHeader)) .* 1000, data(lineNum).dataEmg(end, :) .* 1000, 'LineWidth', 2);  % source data is ~-.010 to +.010   plot in ms and mV

    ax = get(gca); ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
    
%     uistack([hvpz1vpz2, hbackscatter, hvvdd, hvstim, hemg]);
    uistack(hvpz1, 'top')
%     uistack(hvpz2, 'top')
%     uistack(hvpz1vpz2, 'top')
    uistack(hbackscatter, 'top')
    uistack(hvvdd, 'top')
    uistack(hvstim, 'top')
    uistack(hemg, 'top')
    
    lgd12 = legend([hvpz1, hbackscatter, hvvdd, hvstim, hemg]);
    legend('piezo (V)', 'backscatter (mV)', 'vdd (V)', 'stim (V)', 'EMG (mV)');
    legend boxoff;
    lgd12.FontSize = 12; lgd12.EdgeColor = [.65, .65, .65];
    lgd12.Location = 'eastoutside';

    xlabel('time (ms)')
    ylabel('voltage')
    title('Wireless stimulation waveforms')
    xlim([-1, 13]); 
    ylim([-4.2, 4.2]);
    figure(12); set(gcf, 'Color', 'White', 'pos', [40, 40, 40+figSize(1), 40+figSize(2)]); box off; drawnow; set(gca, 'FontSize', figFontSize);
%     figure(12); saveas(gcf, sprintf('%s\time_backscatter_lineNum_%u.png', figuresFoldername, lineNum))
    
    xlim([-.2, .3]); 
    ylim([-4.2, 4.2]);
    xlim([-.2, (data(lineNum).tdc + data(lineNum).duration + .5)]);
    figure(12); set(gcf, 'Color', 'White', 'pos', [40, 40, 40+figSize(1), 40+figSize(2)]); box off; drawnow; set(gca, 'FontSize', figFontSize);
%     figure(12); saveas(gcf, sprintf('%s\time_backscatter_zoom_lineNum_%u.png', figuresFoldername, lineNum))    


    
    %======== figure for linked time backscatter
    figure(14); hold on
    ax1 = subplot(2,1,1); hold on
    hvstim_2 = plot(data(lineNum).vstimTime .* 1000, data(lineNum).vstim, 'LineWidth', 2);  % source data is ~0 to 1 V   plot in ms and V
    hvvdd_2 = plot(data(lineNum).vstimTime .* 1000, data(lineNum).vvdd, 'LineWidth', 2);  % source data is ~0 to 3 V   plot in ms and V
%     hvpz1vpz2_2 = plot(data(lineNum).vstimTime .* 1000, (data(lineNum).vpz1 - data(lineNum).vpz2));  % time data in s, plot in ms; volt data in volt, plot in volt; plot in ms and V
    hemg_2 = plot((data(lineNum).emgTimeTimeDecimated  - (leftShiftTimeOfFlight + leftShiftHeader)) .* 1000, data(lineNum).dataEmg(end, :) .* 1000, 'LineWidth', 2);  % source data is ~-.010 to +.010   plot in ms and mV
    hvpz1_2 = plot(data(lineNum).vstimTime .* 1000, data(lineNum).vpz1, 'LineWidth', .5);  % source data is ~-3 to 3 V   plot in ms and V
    ylim([-1, 4])
%     set(gca,'XColor','white')
%     set(gca,'xtick',[])
%     set(gca,'xticklabel',[])    
    uistack(hvpz1_2, 'top'); 
    uistack(hvvdd_2, 'top'); 
    uistack(hvstim_2, 'top'); 
    uistack(hemg_2, 'top')
    lgd14_1 = legend('Pz+', 'Vdd2p5', 'Velec', 'Vemg');
    legend 'boxoff';
    lgd14_1.Location = 'eastoutside';
    title('Stimulation waveforms');
    ylabel('Voltage (V)');
    xlabel('Time (ms)');



    ax2 = subplot(2,1,2); hold on
    hbackscatter_2 = plot(backscatterToPlotTime .* 1000, backscatterToPlotVolt .* 1000 );  % source data is limited to ~-2048 to +2048; typically ~-800 to +800; plot in ms and mV
    ylim([-40, 40])
    lgd14_2 = legend('Vbackscatter');
    legend 'boxoff';
    lgd14_2.Location = 'eastoutside';
    ylabel('Voltage (mV)');
    xlabel('Time (ms)');


    
    set(hvpz1_2, 'color', color5{5});  % color5 goes from dark to light
    set(hbackscatter_2, 'color', color5{4});
    set(hvvdd_2, 'color', color5{3});
    set(hvstim_2, 'color', color5{2});
    set(hemg_2, 'color', color5{1});
      
%     set(gcf, 'Color', 'white');
    linkaxes([ax1, ax2],'x')    
    xlim([-.4, .5]); 

    figSize = [900 1100]; % [width height]
    monitorPos = get(groot, 'MonitorPositions');
    monitorPixels = (monitorPos(:,4) - monitorPos(:,2)).*(monitorPos(:,3) - monitorPos(:,1));
    [~,monitorToUse] = max(monitorPixels);  % use largest monitor
    figPos = [monitorPos(monitorToUse, 1) + 40, monitorPos(monitorToUse, 4) - figSize(2) - 60, figSize(1), figSize(2)];
    set(14, 'pos', figPos); 
    
    
    xlim(ax1, [-1, 9.5]); 
    ylim(ax1, [-3.2, 4.0]);
    ylim(ax2, [-130, 130])
    yticks([-40 0 40])
    % save zoom out
    fn_format_and_save_figure(14, [pwd filesep figuresFoldername filesep 'waveforms_zoomOut'], [NaN])

    set(hvpz1_2,'XData', get(hvpz1_2, 'XData') .* 1000);  % rescale X so in us not ms
    set(hbackscatter_2,'XData', get(hbackscatter_2, 'XData') .* 1000);  % rescale X so in us not ms
    set(hvvdd_2,'XData', get(hvvdd_2, 'XData') .* 1000);  % rescale X so in us not ms
    set(hvstim_2,'XData', get(hvstim_2, 'XData') .* 1000);  % rescale X so in us not ms
    set(hemg_2,'XData', get(hemg_2, 'XData') .* 1000);  % rescale X so in us not ms
    xlabel(ax1, 'Time (\mus)')
    xlabel(ax2, 'Time (\mus)')
    xlim(ax1, 1000.*[-.2, .35]); 
    ylim(ax1, [-.8, 3.6]);
    ylim(ax2, [-120, 120])
    % save zoom in
    fn_format_and_save_figure(14, [pwd filesep figuresFoldername filesep 'waveforms_zoomIn'], [NaN])

    delete(hemg_2)
    hvpz1_2.LineWidth = 1.75;
    xlim(ax1, 1000.*[-.013, .010]); 
    ylim(ax1, [-.8, 3.6]);
    ylim(ax2, [-120, 120])
    % save zoom header
    fn_format_and_save_figure(14, [pwd filesep figuresFoldername filesep 'waveforms_zoomHeader'], [NaN])
    
    hbackscatter_2shifted = copyobj(hbackscatter_2,ax2); %copy children to new parent axes i.e. the subplot axes
    set(hbackscatter_2shifted,'XData', get(hbackscatter_2shifted, 'XData') - (data(lineNum).tdc + data(lineNum).usDuration).*1e6 );  % rescale X so in us not ms
    set(hbackscatter_2shifted, 'color', [.3, .3, .3])
    xlim(ax2, [-90, -81.5]); 
    ylim(ax2,[-47, 47])
    hbackscatter_2.LineWidth = 1.15;
    hbackscatter_2shifted.LineWidth = 1.15;
    % save zoom backscatter
    drawnow
    fn_format_and_save_figure(14, [pwd filesep figuresFoldername filesep 'waveforms_zoomBackscatter'], [NaN])

    clear hvstim_2 hvvdd_2 hvpz1_2 hbackscatter_2 hemg_2
    

    
    
    
    
        %======== figure for linked time backscatter
    figure(16); clf; hold on

    figSize = [720 1100]; % [width height]
    monitorPos = get(groot, 'MonitorPositions');
    monitorPixels = (monitorPos(:,4) - monitorPos(:,2)).*(monitorPos(:,3) - monitorPos(:,1));
    [~,monitorToUse] = max(monitorPixels);  % use largest monitor
    figPos = [monitorPos(monitorToUse, 1) + 40, monitorPos(monitorToUse, 4) - figSize(2) - 60, figSize(1), figSize(2)];
    set(16, 'pos', figPos); 
    

    ax1 = subplot(3,1,1); hold on
    ax2 = subplot(3,1,2); hold on
    ax3 = subplot(3,1,3); hold on

    
    
    ax1PosOrig = get(ax1, 'pos');
    ax2PosOrig = get(ax2, 'pos');
    ax3PosOrig = get(ax3, 'pos');
    widthToUse = max([ax1PosOrig(3), ax2PosOrig(3), ax3PosOrig(3)]);
    heightNominal = ax1PosOrig(4);
%     set(ax1, 'pos', [ax1PosOrig(1) (ax1PosOrig(2)+.08) widthToUse .45.*heightNominal]);
%     set(ax2, 'pos', [ax2PosOrig(1) (ax2PosOrig(2)+.16) widthToUse .45.*heightNominal]);
%     set(ax3, 'pos', [ax3PosOrig(1) ax3PosOrig(2) widthToUse 1.2.*heightNominal]);    
    set(ax1, 'pos', [ax1PosOrig(1) (ax1PosOrig(2)+.00) widthToUse 1.*heightNominal]);
    set(ax2, 'pos', [ax2PosOrig(1) (ax2PosOrig(2)+.0) widthToUse 1.*heightNominal]);
    set(ax3, 'pos', [ax3PosOrig(1) ax3PosOrig(2) widthToUse 1.*heightNominal]);  

    set(ax1, 'pos', [ax1PosOrig(1) (ax1PosOrig(2)+.00) widthToUse 1.15.*heightNominal]);
    set(ax2, 'pos', [ax2PosOrig(1) (ax2PosOrig(2)+.0) widthToUse 0.65.*heightNominal]);
    set(ax3, 'pos', [ax3PosOrig(1) ax3PosOrig(2) widthToUse 0.90.*heightNominal]);      
    
    
    axes(ax1)
    hvstim_2 = plot(data(lineNum).vstimTime .* 1000, data(lineNum).vstim, 'LineWidth', 2, 'Parent', ax1);  % source data is ~0 to 1 V   plot in ms and V
    hvvdd_2 = plot(data(lineNum).vstimTime .* 1000, data(lineNum).vvdd, 'LineWidth', 2, 'Parent', ax1);  % source data is ~0 to 3 V   plot in ms and V
%     hvpz1vpz2_2 = plot(data(lineNum).vstimTime .* 1000, (data(lineNum).vpz1 - data(lineNum).vpz2));  % time data in s, plot in ms; volt data in volt, plot in volt; plot in ms and V
    hvpz1_2 = plot(data(lineNum).vstimTime .* 1000, data(lineNum).vpz1, 'LineWidth', .5, 'Parent', ax1);  % source data is ~-3 to 3 V   plot in ms and V
    ylim([-1, 4])  
    uistack(hvpz1_2, 'top'); 
    uistack(hvvdd_2, 'top'); 
    uistack(hvstim_2, 'top'); 
    lgd14_1 = legend('Pz+', 'Vdd2p5', 'Velec');
    legend 'boxoff';
    lgd14_1.Location = 'northeast';
    title('Stimulation waveforms');
    ylabel('Voltage (V)');
    xlabel('Time (ms)');



    axes(ax2)
    hbackscatter_2 = plot(backscatterToPlotTime .* 1000, backscatterToPlotVolt .* 1000, 'LineWidth', .5, 'Parent', ax2);  % source data is limited to ~-2048 to +2048; typically ~-800 to +800; plot in ms and mV
    ylim([-40, 40])
    lgd14_2 = legend('Vbackscatter');
    legend 'boxoff';
    lgd14_2.Location = 'northeast';
    ylabel('Voltage (mV)');
    xlabel('Time (ms)');
    
    
    
    axes(ax3)
    hemg_2 = plot((data(lineNum).emgTimeTimeDecimated  - (leftShiftTimeOfFlight + leftShiftHeader)) .* 1000, data(lineNum).dataEmg(end, :) .* 1000, 'LineWidth', 2, 'Parent', ax3);  % source data is ~-.010 to +.010   plot in ms and mV
    ylim([-40, 40])
    lgd14_3 = legend('Vemg');
    legend 'boxoff';
    lgd14_3.Location = 'northeast';
    ylabel('Voltage (mV)');
    xlabel('Time (ms)');


    
    set(hvpz1_2, 'color', color5{5});  % color5 goes from dark to light
    set(hbackscatter_2, 'color', color5{4});
    set(hvvdd_2, 'color', color5{3});
    set(hvstim_2, 'color', color5{2});
    set(hemg_2, 'color', color5{1});
      
%     set(gcf, 'Color', 'white');
    linkaxes([ax1, ax2, ax3],'x')    
    xlim([-.4, .5]); 


    
    
    xlim(ax1, [-.7, 9.25]); 
    ylim(ax1, [-.2, 3.5]);
    ylim(ax2, [-45, 45]);
    ylim(ax3, [-3.5, 4.5]);
    axes(ax2)
    yticks([-40 0 40])
%     set(ax2, 'YTicksBetween', 4)
    yTick = ax2.YRuler.TickValues;
    yTickSpacing = yTick(end) - yTick(end-1);
    yMinorTickSpacing = yTickSpacing ./ 4;
    yMinorTickValues = (yTick(1) - yTickSpacing):yMinorTickSpacing:(yTick(end) + yTickSpacing);
    yMinorTickValues = yMinorTickValues((yMinorTickValues >= ax2.YRuler.Limits(1)) & (yMinorTickValues < ax2.YRuler.Limits(2)));
    ax2.YRuler.MinorTickValues = yMinorTickValues;



    p1 = plot(xlim(ax1), [0 0], 'color', [.7 .7 .7], 'LineWidth', .15, 'Parent', ax1);
    uistack(p1, 'bottom'); 
    p2 = plot(xlim(ax2), [0 0], 'color', [.7 .7 .7], 'LineWidth', .15, 'Parent', ax2);
    uistack(p2, 'bottom'); 
    p3 = plot(xlim(ax3), [0 0], 'color', [.7 .7 .7], 'LineWidth', .15, 'Parent', ax3);
    uistack(p3, 'bottom'); 
    
%     yticks([-40 0 40])
    % save zoom out
    fn_format_and_save_figure(16, [pwd filesep figuresFoldername filesep 'waveforms_zoomOut'], [NaN])

    set(hvpz1_2,'XData', get(hvpz1_2, 'XData') .* 1000);  % rescale X so in us not ms
    set(hbackscatter_2,'XData', get(hbackscatter_2, 'XData') .* 1000);  % rescale X so in us not ms
    set(hvvdd_2,'XData', get(hvvdd_2, 'XData') .* 1000);  % rescale X so in us not ms
    set(hvstim_2,'XData', get(hvstim_2, 'XData') .* 1000);  % rescale X so in us not ms
    set(hemg_2,'XData', get(hemg_2, 'XData') .* 1000);  % rescale X so in us not ms
    set(p1,'XData', get(p1, 'XData') .* 1000);  % rescale X so in us not ms
    set(p2,'XData', get(p2, 'XData') .* 1000);  % rescale X so in us not ms
    set(p3,'XData', get(p3, 'XData') .* 1000);  % rescale X so in us not ms
    xlabel(ax1, 'Time (\mus)')
    xlabel(ax2, 'Time (\mus)')
    xlabel(ax3, 'Time (\mus)')    

    
    
    set(ax1, 'pos', ax1PosOrig);
    set(ax2, 'pos', ax2PosOrig);
    set(ax3, 'pos', ax3PosOrig);  
    
    ax1Pos = get(ax1, 'pos');
    ax2Pos = get(ax2, 'pos');
    ax3Pos = get(ax3, 'pos');
    widthToUse = max([ax1Pos(3), ax2Pos(3), ax3Pos(3)]);
%     set(ax1, 'pos', [ax1PosOrig(1) (ax1PosOrig(2)-.03) widthToUse 1.225.*heightNominal]);
%     set(ax2, 'pos', [ax2PosOrig(1) (ax2PosOrig(2)-.06) widthToUse .75.*heightNominal]);
%     set(ax3, 'pos', [ax3PosOrig(1) ax3PosOrig(2) widthToUse .45.*heightNominal]);  
    set(ax1, 'pos', [ax1PosOrig(1) (ax1PosOrig(2)+.00) widthToUse 1.*heightNominal]);
    set(ax2, 'pos', [ax2PosOrig(1) (ax2PosOrig(2)+.0) widthToUse 1.*heightNominal]);
    set(ax3, 'pos', [ax3PosOrig(1) ax3PosOrig(2) widthToUse 1.*heightNominal]);  
    set(ax1, 'pos', [ax1PosOrig(1) (ax1PosOrig(2)+.00) widthToUse 1.1.*heightNominal]);
    set(ax2, 'pos', [ax2PosOrig(1) (ax2PosOrig(2)+.0) widthToUse 0.65.*heightNominal]);
    set(ax3, 'pos', [ax3PosOrig(1) ax3PosOrig(2) widthToUse 0.85.*heightNominal]);      
    
    xlim(ax1, 1000.*[-.2, .35]); 
    ylim(ax1, [-.2, 3.5]);
    ylim(ax2, [-45, 45])
    ylim(ax3, [-.65, .1])
    axes(ax3)
    yticks([-.5 0])
    % save zoom in
    fn_format_and_save_figure(16, [pwd filesep figuresFoldername filesep 'waveforms_zoomIn'], [NaN])

    delete(hemg_2)
    hvpz1_2.LineWidth = 1.75;
    xlim(ax1, 1000.*[-.013, .010]); 
    ylim(ax1, [-.8, 3.6]);
    ylim(ax2, [-120, 120])
    % save zoom header
    fn_format_and_save_figure(16, [pwd filesep figuresFoldername filesep 'waveforms_zoomHeader'], [NaN])
    
    hbackscatter_2shifted = copyobj(hbackscatter_2,ax2); %copy children to new parent axes i.e. the subplot axes
    set(hbackscatter_2shifted,'XData', get(hbackscatter_2shifted, 'XData') - (data(lineNum).tdc + data(lineNum).usDuration).*1e6 );  % rescale X so in us not ms
    set(hbackscatter_2shifted, 'color', [.3, .3, .3])
    xlim(ax2, [-90, -81.5]); 
    ylim(ax2,[-47, 47])
    hbackscatter_2.LineWidth = 1.15;
    hbackscatter_2shifted.LineWidth = 1.15;
    % save zoom backscatter
    drawnow
    fn_format_and_save_figure(16, [pwd filesep figuresFoldername filesep 'waveforms_zoomBackscatter'], [NaN])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %======== first histogram figure for comparing recharge pop to stim pop
%     figure(11)  
%     subplot(1,2,1); hold on
%     histogram(data(lineNum).backscatterL1normRecharge)
%     histogram(data(lineNum).backscatterL1normStim)
%     xlabel('backscatter L1 norm')
%     ylabel('occurances')
%     title('distribution of backscatter L1 Norm, post-stim vs. post-recharge')
%     subplot(1,2,2)
%     histogram(data(lineNum).backscatterL1normRatio); hold on
%     histogram(data(lineNum).backscatterL1normRecharge(2:end)./data(lineNum).backscatterL1normRecharge(1:(end-1)))
%     xlabel('backscatter L1 norm ratio')
%     ylabel('occurances')
%     title('distribution of ratio of stim to recharge backscatter L1 norm')
%     
%     figure(11); set(gcf, 'Color', 'White', 'pos', [40, 40, 40+figSize(1), 40+figSize(2)]); box off; drawnow; %set(gca, 'FontSize', figFontSize);
% %     figure(11); saveas(gcf, sprintf('%s\backscatter_comparison_lineNum_%u.png', figuresFoldername, lineNum))    

    bsRR = data(lineNum).backscatterL1normRecharge(2:end)./data(lineNum).backscatterL1normRecharge(1:(end-1));
    bsRS = data(lineNum).backscatterL1normRatio(1:(end-1));

    
    
    
    %======== second histogram figure ratiometric regions
    figure(13); hold on
    movingAverageWindow = 1;
    bsRR = conv(bsRR, (ones(1,movingAverageWindow)./movingAverageWindow), 'same');
    bsRR = bsRR((movingAverageWindow+0):(end - movingAverageWindow));
    bsRS = conv(bsRS, (ones(1,movingAverageWindow)./movingAverageWindow), 'same');
    bsRS = bsRS((movingAverageWindow+0):(end - movingAverageWindow));
    
    histogram(bsRR, 'BinWidth', .003)
    histogram(bsRS, 'BinWidth', .003)
    xLimUse = xlim();
    
    distributionX = ((xLimUse(1) - 1.4):.0001:(xLimUse(2) + 1.4));
    distributionRRadd = normpdf(distributionX, mean(bsRR), std(bsRR));
    distributionRSadd = normpdf(distributionX, mean(bsRS), std(bsRS));  %length(bsRR).*
%     distributionRR = normpdf(distributionX, mean(bsRR), std(bsRR));
%     distributionRS = normpdf(distributionX, mean(bsRS), std(bsRS));
%     distributionRR = normpdf(distributionX, 0, 1);
%     distributionRS = normpdf(distributionX, 0, 1);
    
    plot(distributionX, distributionRRadd)
    plot(distributionX, distributionRSadd)
%     xlim(xLimUse)


    %======== third histogram figure additive comparison ratiometric
    %regions recharge minus stim all over recharge
    figure(15); clf; hold on
    beginningSkip = 1;
    endSkip = 1;
    % next recharge minus this recharge all over this recharge
    bsRRadd = (data(lineNum).backscatterL1normRecharge(beginningSkip+2:end-endSkip) - data(lineNum).backscatterL1normRecharge(beginningSkip+1:(end-endSkip-1))) ./ data(lineNum).backscatterL1normRecharge(beginningSkip+1:(end-endSkip-1));
    % this stim minus this recharge all over this recharge
    bsRSadd = (data(lineNum).backscatterL1normStim(beginningSkip+1:(end-endSkip)) - data(lineNum).backscatterL1normRecharge(beginningSkip+1:(end-endSkip))) ./ data(lineNum).backscatterL1normRecharge(beginningSkip+1:(end-endSkip));
    
    
    
    movingAverageWindow = 1;
    bsRRadd = conv(bsRRadd, (ones(1,movingAverageWindow)./movingAverageWindow), 'same');
    bsRRadd = bsRRadd((movingAverageWindow+0):(end - movingAverageWindow));
    bsRSadd = conv(bsRSadd, (ones(1,movingAverageWindow)./movingAverageWindow), 'same');
    bsRSadd = bsRSadd((movingAverageWindow+0):(end - movingAverageWindow));
    
    histogram(bsRRadd, 'BinWidth', .00195, 'EdgeColor', [0.15 0.15 0.15], 'FaceColor', color5{2}, 'FaceAlpha', 1.0, 'LineWidth', 0.3)
    histogram(bsRSadd, 'BinWidth', .00195, 'EdgeColor', [0.15 0.15 0.15], 'FaceColor', color5{end-1}, 'FaceAlpha', 1.0, 'LineWidth', 0.3)
    xLimUse = xlim();
    
    distributionX = ((xLimUse(1) - 1.4):.0001:(xLimUse(2) + 1.4));
    distributionRRadd = normpdf(distributionX, mean(bsRRadd), std(bsRRadd));
    distributionRSadd = normpdf(distributionX, mean(bsRSadd), std(bsRSadd));  %length(bsRRadd).*
%     distributionRR = normpdf(distributionX, mean(bsRRadd), std(bsRRadd));
%     distributionRS = normpdf(distributionX, mean(bsRSadd), std(bsRSadd));
%     distributionRR = normpdf(distributionX, 0, 1);
%     distributionRS = normpdf(distributionX, 0, 1);
    
%     plot(distributionX, distributionRRadd)
%     plot(distributionX, distributionRSadd)
    xlim([-.11, 0.05])
    title('Through-tissue backscatter state differentiation')
    legend(' device simulated non-operating',' device operating')
    legend boxoff
    xlabel('backscatter relative difference')
    ylabel(sprintf('count out of N = %u', length(bsRRadd)))
    fn_format_and_save_figure(15, [pwd filesep figuresFoldername filesep 'histogram_L2_1mov'], [900*.75 550*.75])
    
    
    % find equal point
%     [M,I] = min(abs(distributionRRadd - distributionRSadd));
%     bModDiscriminateThreshold = distributionX(I)
    bModDiscriminateThreshold = -.0451;
    bModDiscriminateThreshold = -.0449995;
    BER_falseZero = normcdf(bModDiscriminateThreshold, mean(bsRRadd), std(bsRRadd))
    BER_falseOne = 1 - normcdf(bModDiscriminateThreshold, mean(bsRSadd), std(bsRSadd))
    BER_average = mean([BER_falseZero, BER_falseOne]);
    BER_average    
    
    bModDiscriminateThreshold = mean([mean(bsRRadd), mean(bsRSadd)]);
    BER = 1 - normcdf(bModDiscriminateThreshold, mean(bsRSadd), std(bsRSadd));
    BER
    
    
    mean(bsRSadd)
    std(bsRSadd)
    mean(bsRRadd)
    std(bsRRadd)
end





