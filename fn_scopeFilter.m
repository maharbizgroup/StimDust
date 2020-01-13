function [backscatterTime, backscatterFilt, backscatterL1normRecharge, backscatterL1normStim, backscatterL1normRatio] = fn_scopeFilter(rawBackscatter, freqSampling, freqCarrier, tdcDuration, usDuration, backscatterWindow)

%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2018-04-04
% All rights reserved.
%========================================


            
            % create vstim and vdd low-pass filter and apply filter
            firNumTaps = 4690;
            cutoffFreq = 1400000;   % Hz
            scopeSamplingFreq = 1./scopeSampleTimeStep;
            nyquistRate = scopeSamplingFreq./2;  
            Wn = [cutoffFreq./nyquistRate];
            firCoefsScopeLowPass = fir1(firNumTaps, Wn, 'low');
%             vstim = filtfilt(firCoefsScopeLowPass, 1, vstim);
%             vvdd = filtfilt(firCoefsScopeLowPass, 1, vvdd);

            % create pz1 pz2 low-pass filter and apply filter
            firNumTaps = 4690;
            cutoffFreq = 4400000;   % Hz
            scopeSamplingFreq = 1./scopeSampleTimeStep;
            nyquistRate = scopeSamplingFreq./2;  
            Wn = [cutoffFreq./nyquistRate];
            firCoefsScopeLowPass = fir1(firNumTaps, Wn, 'low');
%             vpz1 = filtfilt(firCoefsScopeLowPass, 1, vpz1);
%             vpz2 = filtfilt(firCoefsScopeLowPass, 1, vpz2);


    PLOTTING = 0;

    persistent firCoefsBackscatterBandPass
    if(isempty(firCoefsBackscatterBandPass))
        %======== backscatter filter
        firNumTaps = 456;
        relFreqWindow = [0.31 1.4];
        signalFreq = freqCarrier;  % Hz
        %lowerCutoffFreq = 200;  % Hz
        %upperCutoffFreq = 10000;   % Hz
        nyquistRate = freqSampling./2;  
        Wn = (signalFreq.*relFreqWindow)./nyquistRate;
        firCoefsBackscatterBandPass = fir1(firNumTaps, Wn, 'bandpass');
        % freqz(firCoefsRawHighPass, 1, 1000);
    end
    

    backscatterTime = (1:length(rawBackscatter)) ./ freqSampling;
    
    backscatterWindow = [30e-6, 50e-6];
    backscatterWindowSamplesRecharge = floor(backscatterWindow .* freqSampling);
    backscatterWindowSamplesStim = backscatterWindowSamplesRecharge + floor((tdcDuration + usDuration) .* freqSampling);
    windowDCOffset = [floor((tdcDuration + 20e-6) .* freqSampling), floor((tdcDuration + usDuration - 20e-6) .* freqSampling)];
    windowDCOffset = [floor(20e-6 .* freqSampling), floor(60e-6 .* freqSampling)];
    

    dataBackscatter = rawBackscatter;
    dataBackscatter = dataBackscatter - mean(dataBackscatter(windowDCOffset(1):windowDCOffset(2)));
    backscatterFilt = filtfilt(firCoefsBackscatterBandPass, 1, dataBackscatter);

    if(PLOTTING)
        figure(10); hold off
        subplot(1,2,1)
        plot(1:length(dataBackscatter), dataBackscatter, 'b-');
%         xlim([0, 40e-6]);
%         ylim([-900, 900]); 

        subplot(1,2,2); hold on
        plot(backscatterTime, backscatterFilt, 'b-');
        plot(backscatterTime(backscatterWindowSamplesRecharge), backscatterFilt(backscatterWindowSamplesRecharge), 'r-')
        plot(backscatterTime(backscatterWindowSamplesStim), backscatterFilt(backscatterWindowSamplesStim), 'r-'); hold off
%         xlim([(tdcDuration + usDuration), (tdcDuration + usDuration) + 40e-6]);
%         ylim([-900, 900]); 
    end

    backscatterL1normRecharge = sum(abs(backscatterFilt(backscatterWindowSamplesRecharge(1):backscatterWindowSamplesRecharge(2))));
    backscatterL1normStim = sum(abs(backscatterFilt(backscatterWindowSamplesStim(1):backscatterWindowSamplesStim(2))));
    backscatterL1normRatio = backscatterL1normStim ./ backscatterL1normRecharge;

end