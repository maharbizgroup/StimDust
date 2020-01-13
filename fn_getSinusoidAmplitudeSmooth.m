function [amplitude, pkpk, rms] = fn_getSinusoidAmplitudeSmooth(inputROI, freqSampling, freqCarrier)
    
%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2018-04-24
% All rights reserved.
%========================================


%======== persistent version
%     persistent firCoefsRawLowPass
%     if(isempty(firCoefsRawLowPass))
%         % raw data low-pass filter
%         firNumTaps = 6.*freqSampling.*freqCarrier;
%         firNumTaps = 800;
%         cutoffFreq = freqCarrier .* 6;   % Hz
%         nyquistRate = freqSampling./2;  
%         Wn = [cutoffFreq./nyquistRate];
%         firCoefsRawLowPass = fir1(firNumTaps, Wn, 'low');
%         % freqz(firCoefsRawLowPass, 1, 1000);
%     end
    

    %======== non persistent version
    % raw data low-pass filter
    firNumTaps = 6.*freqSampling.*freqCarrier;
    firNumTaps = 800;
    if length(inputROI) < 3.*firNumTaps
        firNumTaps = floor(length(inputROI) ./ 3 - 1);
    end
    cutoffFreq = freqCarrier .* 6;   % Hz
    nyquistRate = freqSampling./2;  
    Wn = [cutoffFreq./nyquistRate];
    firCoefsRawLowPass = fir1(firNumTaps, Wn, 'low');
    % freqz(firCoefsRawLowPass, 1, 1000);
    
    y_lowpass = filtfilt(firCoefsRawLowPass, 1, inputROI);

    amplitude = max(y_lowpass);
    pkpk = max(y_lowpass) - min(y_lowpass);
    rms = sqrt(mean((y_lowpass - mean(y_lowpass)).^2));
    
%     figure(251); hold on
%     plot(inputROI)
%     plot(y_lowpass)
end
