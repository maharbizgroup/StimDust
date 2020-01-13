function [backscatterTime, backscatterFilt, backscatterL1normRecharge, backscatterL1normStim, backscatterL1normRatio] = fn_stimBackscatterProcess(rawBackscatter, freqSampling, freqCarrier, tdcDuration, usDuration, backscatterWindow, plotOn)

%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2019-02-01
% All rights reserved.
%========================================


    PLOTTING = plotOn;
    PLOTTING = 0;
    FFTON = 0;

    persistent timesCalled
    if isempty(timesCalled)
        timesCalled = 0;
    end
    timesCalled = timesCalled + 1;
    if(mod(timesCalled, 10) == 0)
        PLOTTING=1;
    else
        PLOTTING=0;
    end
    
    persistent firCoefsBackscatterBandPass
    if(isempty(firCoefsBackscatterBandPass))
        %======== backscatter filter
        firNumTaps = 456;
        relFreqWindow = [0.95 1.08];
        signalFreq = freqCarrier;  % Hz
        %lowerCutoffFreq = 200;  % Hz
        %upperCutoffFreq = 10000;   % Hz
        nyquistRate = freqSampling./2;  
        Wn = (signalFreq.*relFreqWindow)./nyquistRate;
        firCoefsBackscatterBandPass = fir1(firNumTaps, Wn, 'bandpass');
        % freqz(firCoefsRawHighPass, 1, 1000);
    end
    

    backscatterTime = (1:length(rawBackscatter)) ./ freqSampling;
    
    backscatterWindowIndicesRecharge = floor(backscatterWindow(1) .* freqSampling):ceil(backscatterWindow(2) .* freqSampling);
    backscatterWindowIndicesStim = backscatterWindowIndicesRecharge + floor((tdcDuration + usDuration) .* freqSampling);
    windowDCOffset = [floor((tdcDuration + 20e-6) .* freqSampling), floor((tdcDuration + usDuration - 20e-6) .* freqSampling)];
    windowDCOffset = [floor(5e-6 .* freqSampling), floor(tdcDuration .* freqSampling)];
    

    dataBackscatter = rawBackscatter;
    if(windowDCOffset(2) <= length(rawBackscatter))
        dataBackscatter = dataBackscatter - mean(dataBackscatter(windowDCOffset(1):windowDCOffset(2)));
    end
    backscatterFilt = filtfilt(firCoefsBackscatterBandPass, 1, dataBackscatter);
    
    if((backscatterWindowIndicesStim(2) > length(rawBackscatter)) || (windowDCOffset(2) > length(rawBackscatter)))
        backscatterL1normRecharge = NaN;
        backscatterL1normStim = NaN;
        backscatterL1normRatio = NaN;
    else

        if(PLOTTING)
            figure(21); hold off; clf; hold on
            %figure; hold off; clf; hold on
%             set(gcf, 'Position', get(0,'Screensize'));
            subplot(2,1,1); hold on
%             plot(1:length(dataBackscatter), dataBackscatter, 'b-');
            plot(backscatterTime, dataBackscatter, 'b-');
            plot(backscatterTime(backscatterWindowIndicesRecharge), dataBackscatter(backscatterWindowIndicesRecharge), 'r-')
            plot(backscatterTime(backscatterWindowIndicesStim), dataBackscatter(backscatterWindowIndicesStim), 'r-'); hold off
            
    %         xlim([0, 40e-6]);
            ylim([-.07, .07]); 

%             figure(22); hold off; clf; hold on
            subplot(2,1,2); hold on
            plot(backscatterTime, backscatterFilt, 'b-');
            plot(backscatterTime(backscatterWindowIndicesRecharge), backscatterFilt(backscatterWindowIndicesRecharge), 'r-')
            plot(backscatterTime(backscatterWindowIndicesStim), backscatterFilt(backscatterWindowIndicesStim), 'r-'); hold off
    %         xlim([(tdcDuration + usDuration), (tdcDuration + usDuration) + 40e-6]);
            ylim([-.1, .1]); 
            ylim([-.03, .03]); 

            if(FFTON)
                hFig = fn_plot_fft(dataBackscatter(backscatterWindowIndicesRecharge), freqSampling);
                figure(21);
                h3 = subplot(1,3,3); hold on
                copyobj(allchild(get(hFig,'CurrentAxes')),h3);
            end
            drawnow
        end

%         backscatterL1normRecharge = sum(abs(backscatterFilt(backscatterWindowIndicesRecharge)));
%         backscatterL1normStim = sum(abs(backscatterFilt(backscatterWindowIndicesStim)));
        norm_p = 2;
        backscatterL1normRecharge = norm(backscatterFilt(backscatterWindowIndicesRecharge), norm_p);
        backscatterL1normStim = norm(backscatterFilt(backscatterWindowIndicesStim), norm_p);
        backscatterL1normRatio = backscatterL1normStim ./ backscatterL1normRecharge;
    end

end