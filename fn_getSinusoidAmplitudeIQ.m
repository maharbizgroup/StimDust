function [amplitude, searchRMS] = fn_getSinusoidAmplitudeIQ(inputROI, freqSampling, freqCarrier, findRegion)
%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2019-06-26
% All rights reserved.
%========================================


% Return amplitude of sinusoidal signal using IQ demodulation
%
% Inputs:
%    inputROI: input signal
%    freqSampling: 
%    freqCarrier: freq of sinusoid
%    findRegion: bool; automatically find relevant region (usually no)
%
% Outputs:
%    amplitude: amplitude of whole inputROI
%    searchRMS: amplitude of region found if findRegion == 1
%
% Author: David Piech
% Date: 2018-03-26
% 

    PLOTTING = 1;

    inputROIorig = inputROI;
    inputROIBoxcar = conv(inputROI, ones(1,4)./4, 'same');
    inputROI = inputROIBoxcar;
    
    decimationFactor = floor((freqSampling./freqCarrier) ./ 30);
    freqSampling = freqSampling ./ decimationFactor;

    inputROI = inputROI - mean(inputROI);
    inputROI = inputROI(1:decimationFactor:end);

    mI = cos(2 * pi * (freqCarrier / (freqSampling)) * (1:length(inputROI)));
    mQ = -sin(2 * pi * (freqCarrier / (freqSampling)) * (1:length(inputROI)));

    yI = inputROI .* mI;
    yQ = inputROI .* mQ;
        
    %======== filter persistent version
%     clear firCoefsRawLowPass
%     persistent firCoefsRawLowPass
%     if(isempty(firCoefsRawLowPass))
%         % raw data low-pass filter
%         firNumTaps = floor(2.*6.*freqSampling./freqCarrier);
% %         firNumTaps = 100;
%         firNumTaps = floor(min([firNumTaps, length(inputROI)./3]));
%         firNumTaps = 100;
%         if length(inputROI) < 3.*firNumTaps
%             firNumTaps = floor(length(inputROI) ./ 3 - 1)
%         end
%         cutoffFreq = freqCarrier ./ 8;   % Hz
%         nyquistRate = freqSampling./2;  
%         Wn = [cutoffFreq./nyquistRate];
%         firCoefsRawLowPass = fir1(firNumTaps, Wn, 'low');
%         % freqz(firCoefsRawLowPass, 1, 1000);
%     end
    
    %======== filter non-persistent version
    % raw data low-pass filter
    firNumTaps = floor(2.*6.*freqSampling./freqCarrier);
%         firNumTaps = 100;
    firNumTaps = floor(min([firNumTaps, length(inputROI)./3]));
    firNumTaps = 100;
    if length(inputROI) < 3.*firNumTaps
        firNumTaps = floor(length(inputROI) ./ 3 - 1);
    end
    cutoffFreq = freqCarrier ./ 8;   % Hz
    nyquistRate = freqSampling./2;  
    Wn = [cutoffFreq./nyquistRate];
    firCoefsRawLowPass = fir1(firNumTaps, Wn, 'low');
    % freqz(firCoefsRawLowPass, 1, 1000);
    
%     fprintf('wtf %d %d\n', length(firCoefsRawLowPass), length(yI))
    yI_lowpass = 2.*filtfilt(firCoefsRawLowPass, 1, yI);
    yQ_lowpass = 2.*filtfilt(firCoefsRawLowPass, 1, yQ);
    
    yA = 1.*(yI_lowpass.^2 + yQ_lowpass.^2).^(.5);

    if ~findRegion
        amplitude = mean(yA);
        searchRMS = sqrt(mean(inputROIBoxcar.^2));
    else
       
        threshold = .7.*max(yA);
        indicesAboveThreshold = find(yA>threshold);
        ROIIndexmin = min(indicesAboveThreshold);
        ROIIndexmax = max(indicesAboveThreshold);
        ROIwidth = ROIIndexmax - ROIIndexmin;
        
        % use .05 right side for large transducer and .70 right side for small
        % transducer
        
        %== for ztrans large
%         ROI_of_pulse_portion = [0.74, 0.9035];
%         %ROI_of_pulse_portion = [0.82*.4, 0.935*.4];
%         ROIyA_indices = ceil(ROIIndexmin + floor(ROIwidth.*ROI_of_pulse_portion(1))):floor(ROIIndexmin + floor(ROIwidth.*ROI_of_pulse_portion(2)));
        
        %== for ztrans small
%         ROI_of_pulse_indices = [790, 890];
%         ROIyA_indices = ceil(ROIIndexmin + ROI_of_pulse_indices(1)):floor(ROIIndexmin + ROI_of_pulse_indices(2));
        
        %== for xtrans small
%         ROIyA_indices = 1750:2250;
        
        %== for xtrans large
        ROIyA_indices = 1750:2000;
  
        
        ROIyA = yA(ROIyA_indices);  % for backscatters that have no standing wave formed
%         amplitude = min(ROIyA);
        amplitude = mean(ROIyA);
%         ROIyA_sort = sort(ROIyA);
%         amplitude = mean(ROIyA(1:floor(.5.*length(ROIyA_sort))));
%         disp(amplitude)
        
        ROIBoxcar = inputROIBoxcar(ceil(ROIIndexmin + ROIwidth.*.05):floor(ROIIndexmax - ROIwidth.*.70));  % for backscatters that have no standing wave formed
%         plot(ceil(ROIIndexmin + ROIwidth.*.05):1:ceil(ROIIndexmin + ROIwidth.*.55)+length(ROIyA)-1, ROIBoxcar)

        
        if(PLOTTING)
            figure(601); clf; hold on
            plot(inputROI, 'g')
            plot(yA, 'b')
            plot(indicesAboveThreshold, yA(indicesAboveThreshold), 'c')
            plot(ROIyA_indices, ROIyA, 'r')
    %         plot(yI)
    %         plot(yQ)
    %         plot(yI_lowpass)
    %         plot(yQ_lowpass)
        end
        
        
        searchRMS = sqrt(mean(ROIBoxcar.^2));
        
    end
end
