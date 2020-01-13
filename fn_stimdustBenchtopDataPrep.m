function [Misalignment, AngleOffset, AmpPz1Out, Pkpk, Ampvdd, Stimval, BackscatterL1normratio, BackscatterL1normRecharge, BackscatterL1normStim, IndicesRealData] = fn_stimdustBenchtopDataPrep(filenames, scopeChannels, SWEEPTYPE, ROIIndices)

%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2019-05-12
% All rights reserved.
%========================================



PLOTTING = 1;

pz1Ch = scopeChannels.pz1;
pz2Ch = scopeChannels.pz2;
vddCh = scopeChannels.vdd;
stimCh = scopeChannels.stim;

clear pulsesToProcess
% pulsesToProcess = [1 160];

for filenum = 1:size(filenames,1)

        fprintf('\n%s',filenames{filenum})
        filenameNew = ['E:\stimDustData\backscatter_and_scope_data\', filenames{filenum}, '.mat'];
        filenameNew = ['D:\yesfh\stimDustData\backscatter_and_scope_data\', filenames{filenum}, '.mat'];
%         if (~exist('filename', 'var') || ~strcmp(filenameNew, filename))
        clear dataStore
        filename = filenameNew;
        fprintf(' loading...')
        load(filename)
        fprintf(' loaded\n')
%         end


        % stim low-pass filter
        firNumTaps = 800;
        cutoffFreq = .5e6;   % Hz
        nyquistRate = (1 ./ (dataStore(1).scope.xData(1,2) - dataStore(1).scope.xData(1,1)))./2;  
        Wn = cutoffFreq./nyquistRate;
        firCoefsRawLowPass = fir1(firNumTaps, Wn, 'low');
        % freqz(firCoefsRawLowPass, 1, 1000);


%         if exist('pulsesToProcess')
%             pulsesToUse = pulsesToProcess(1):pulsesToProcess(2);
%         else
% %             pulsesToUse = 1:(length(dataStore) - 1);
%             pulsesToUse = 10:(length(dataStore) - 20);
%         end


        %== FOR MISALIGNMENT
%         pulsesToUse = 20:(length(dataStore) - 20);
        pulsesToUse = 1:(length(dataStore)-1);

        
        for m = pulsesToUse
            if(mod(m,10) == 0); fprintf(' %u',m); end

            ROI = dataStore(m).scope.yData(2,ROIIndices);
            if isfield(dataStore(m), 'stepperCurrOffsetMM')
                misalignment(m) = dataStore(m).stepperCurrOffsetMM;
            else
                misalignment(m) = 0;
            end
            pkpk(m) = max(ROI) - min(ROI);
            freqSampling = 1 ./ (dataStore(m).scope.xData(1,2) - dataStore(m).scope.xData(1,1));
            freqCarrier = 1.8e6;  % Hz

            % get amplitudes for each pulse
            % pz plus amp for hydrophone
            
            [ampPz1IQSearch(m), searchRMS(m)] = fn_getSinusoidAmplitudeIQ(dataStore(m).scope.yData(pz1Ch,1:end), freqSampling, freqCarrier, true);
            
            [ampPz1IQ(m), searchRMS(m)] = fn_getSinusoidAmplitudeIQ(dataStore(m).scope.yData(pz1Ch,ROIIndices), freqSampling, freqCarrier, false);
            [ampPz1Smooth(m), pkpkPz1Smooth(m), rmsPz1(m)] = fn_getSinusoidAmplitudeSmooth(dataStore(m).scope.yData(pz1Ch,ROIIndices), freqSampling, freqCarrier);
            if(PLOTTING)
                figure(253); clf; hold on; 
                plot(1:length(dataStore(m).scope.yData(pz1Ch,:)), dataStore(m).scope.yData(pz1Ch,:), 'b-')
                plot(ROIIndices, dataStore(m).scope.yData(pz1Ch,ROIIndices), 'r-')
            end

            [ampPz2IQ(m), searchRMS(m)] = fn_getSinusoidAmplitudeIQ(dataStore(m).scope.yData(pz2Ch,ROIIndices), freqSampling, freqCarrier, false);
            [ampPz2Smooth(m), pkpkPz2Smooth(m), rmsPz2] = fn_getSinusoidAmplitudeSmooth(dataStore(m).scope.yData(pz2Ch,ROIIndices), freqSampling, freqCarrier);        

            [ampPz1Pz2IQ(m), searchRMS(m)] = fn_getSinusoidAmplitudeIQ((dataStore(m).scope.yData(pz1Ch,ROIIndices) - dataStore(m).scope.yData(pz2Ch,ROIIndices)), freqSampling, freqCarrier, false);
            [ampPz1Pz2Smooth(m), pkpkPz1Pz2Smooth(m), rmsPz1Pz2] = fn_getSinusoidAmplitudeSmooth((dataStore(m).scope.yData(pz1Ch,ROIIndices) - dataStore(m).scope.yData(pz2Ch,ROIIndices)), freqSampling, freqCarrier);        

            ampVdd(m) = mean(dataStore(m).scope.yData(vddCh,ROIIndices));

            % get backscatter for each pulse
            backscatterScale = (0.03./0.0896).*0.8/4096.0;   % 0.8 V range divided into 2^12 steps with front-end gain of .0896/.03
            [~, ~, backscatterL1normRecharge(m), backscatterL1normStim(m), backscatterL1normRatio(m)] = fn_stimBackscatterProcess(dataStore(m).backscatter.rawSignal .* backscatterScale, dataStore(m).backscatter.samplingFreq, 1.85e6, 100e-6, 90e-6, [10e-6, 16e-6], 1);
            backscatterBit(m) = abs(backscatterL1normRatio(m) - 1) > .2;

%             fprintf('%d  |  %s\n', m, dataStore(m).backscatter.timeString)
            
            % get ground-truth stim yes/no for each pulse
            stim = dataStore(m).scope.yData(stimCh,:);

            stim_lowpass = filtfilt(firCoefsRawLowPass, 1, stim);

            stim_lowpass(stim_lowpass < .05) = 0;
            stimVal(m) = sum(stim_lowpass) .* (1./freqSampling);  %V*s
            stimBit(m) = stimVal(m) > 4e-5;  % set threshold here

            backscatterCorrect(m) = backscatterBit(m) == stimBit(m);

            if(PLOTTING)
                figure(3)
                plot((dataStore(m).scope.yData(4,:))); hold off
            end

        end

        % align misalignment to peak
        if strcmp(SWEEPTYPE, 'trans') || strcmp(SWEEPTYPE, 'angle')
            useForPeak = pkpkPz1Pz2Smooth;
            [maxPeak, indexPeak] = max(useForPeak);
            misalignment = misalignment - misalignment(indexPeak(1));
        end

        if strcmp(SWEEPTYPE, 'angle')
            if filenum == 2
                misalignment = misalignment - -1.5;
            end
            if filenum == 3
                misalignment = misalignment - -0;
            end

            indicesToUse = find(abs(misalignment) < 2.75);
            k = find(misalignment(indicesToUse) < -2.65);
            if length(k) > 1
                indicesToUse(k(2:end)) = [];
            end
            k = find(misalignment(indicesToUse) > 2.65);
            if length(k) > 1
                indicesToUse(k(2:end)) = [];
            end
        elseif strcmp(SWEEPTYPE, 'ztrans')
            indicesToUse = 1:length(misalignment);
            indicesToUse = 11:3:length(misalignment);
            indicesToUse(pkpkPz1Pz2Smooth(indicesToUse) < .5) = indicesToUse(pkpkPz1Pz2Smooth(indicesToUse) < .5) + 1;
        elseif strcmp(SWEEPTYPE, 'xtransHydrophone')
            indicesToUse = 2:2:length(misalignment);
        else
%             indicesToUse = 1:length(misalignment);
            indicesToUse = 2;
        end
        
        
        
        
        Misalignment(filenum, 1:length(indicesToUse)) = misalignment(indicesToUse);
        AngleOffset(filenum, 1:length(indicesToUse)) = filenames{filenum,2} .* ones(1,length(misalignment(indicesToUse)));
        PzPlusAmp(filenum, 1:length(indicesToUse)) = ampPz1IQ(indicesToUse);
        Pkpk(filenum, 1:length(indicesToUse)) = pkpkPz1Pz2Smooth(indicesToUse);
        Ampvdd(filenum, 1:length(indicesToUse)) = ampVdd(indicesToUse);
        Stimval(filenum, 1:length(indicesToUse)) = stimVal(indicesToUse);
        BackscatterL1normratio(filenum, 1:length(indicesToUse)) = backscatterL1normRatio(indicesToUse);
        BackscatterL1normRecharge(filenum, 1:length(indicesToUse)) = backscatterL1normRecharge(indicesToUse);
        BackscatterL1normStim(filenum, 1:length(indicesToUse)) = backscatterL1normStim(indicesToUse);
        IndicesRealData(filenum, :) = [1, length(indicesToUse)];
%         AmpPz1Out(filenum, 1:length(indicesToUse)) = ampPz1IQ(indicesToUse);
        AmpPz1Out(filenum, 1:length(indicesToUse)) = ampPz1IQSearch(indicesToUse);


        if strcmp(SWEEPTYPE, 'ztransHydrophoneLarge') | strcmp(SWEEPTYPE, 'ztransHydrophoneSmall')
            Misalignment = misalignment;
            AmpPz1Out = ampPz1IQSearch;
        else
        
        %     
    %     
    % %         figure(1); hold off
    % %         plot(dataStore(m).scope.xData(2,:), dataStore(m).scope.yData(2,:))
    % %         ylim([-1 4]); 
    % %         drawnow
    % %         figure(2); hold on
    % %         figure
    % %         plot(misalignment, pkpk); 
    % %         ylim([0 5]); 
    % %         drawnow
        if(PLOTTING)
            figure(250); clf; hold on
            plot(misalignment(indicesToUse), ampPz1IQSearch(indicesToUse), 'b.-', 'DisplayName', 'ampPz1IQ')
%             plot(misalignment(indicesToUse), ampPz1Smooth(indicesToUse), 'DisplayName', 'ampPz1Smooth')
%             plot(misalignment(indicesToUse), pkpkPz1Smooth(indicesToUse), 'DisplayName', 'pkpkPz1Smooth')
%             plot(misalignment(indicesToUse), ampPz2IQ(indicesToUse), 'DisplayName', 'ampPz2IQ')
%             plot(misalignment(indicesToUse), ampPz2Smooth(indicesToUse), 'DisplayName', 'ampPz2Smooth')
%             plot(misalignment(indicesToUse), pkpkPz2Smooth(indicesToUse), 'DisplayName', 'pkpkPz2Smooth')
% 
%             plot(misalignment(indicesToUse), ampPz1Pz2IQ(indicesToUse), 'DisplayName', 'ampPz1Pz2IQ')
%             plot(misalignment(indicesToUse), ampPz1Pz2Smooth(indicesToUse), 'DisplayName', 'ampPz1Pz2Smooth')
%             plot(misalignment(indicesToUse), pkpkPz1Pz2Smooth(indicesToUse), 'DisplayName', 'pkpkPz1Pz2Smooth')
%             plot(misalignment(indicesToUse), ampVdd(indicesToUse), 'DisplayName', 'ampVdd')
% 
% 
%             plot(misalignment(indicesToUse), backscatterL1normRatio(indicesToUse), 'DisplayName', 'backscatterL1normRatio')
%             plot(misalignment(indicesToUse), stimVal(indicesToUse) .* 1e4, 'DisplayName', 'stimVal')
            legend
        %         figure
        %         plot(1:length(backscatterL1normRatio), backscatterL1normRatio)
        end

        
        if (filenum ~= size(filenames,1))
            clear misalignment pkpk ampPz1IQ ampPz1Smooth pkpkPz1Smooth ampPz2IQ ampPz2Smooth pkpkPz2Smooth ampPz1Pz2IQ ampPz1Pz2Smooth pkpkPz1Pz2Smooth ampVdd backscatterL1normRecharge backscatterL1normStim backscatterL1normRatio backscatterBit stimVal stimBit backscatterCorrect ampPz1IQSearch
        end
    end

end