function [intensity_timeAverage, p_amplitude, MI] = fn_acousticMeasurementsHydrophone(vAmpHydrophone, distance)
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

% Return acoustic intensity, pressure from hydrophone voltage
%
% Inputs:
%    vAmpHydrophone: amplitude of raw hydrophone recording (note: not pk-pk nor rms)
%    distance: distance from ; only used when computing 'derated' values
%
% Outputs:
%    intensity_timeAverage: acoustic intensity in W/m^2, averaged over the whole input
%    p_amplitude: pressure amplitude in Pa
%    MI: mechanical index
%
% Internal switches:
%    USE_DERATED derates the acoustic intensity according to the distance and tissue attenuation
%    derated_dBpcmpMHz: the tissue attenuation
%
% Author: David Piech
% Date: 2018-04-18
% 
USE_DERATED = 1;
PLOTTING = 1;

freq_ultrasound = 1.85e6;  % Hz
hydrophone_sensitivity = 251.189e-9;  % V/Pa
amplifier_amplitude_ratio = 10; 
density_water = 1000;  % kg/m^3
speed_sound_water = 1540;  % m/s

vRms = vAmpHydrophone .* 1/sqrt(2);

% derated_dBpcmpMHz = 0;
% derated_dBpcmpMHz = 2.2 / 1.85;  % for muscle at 1.8 MHz
% derated_dBpcmpMHz = 2.2 / 1.8;  % for castor oil at 1.8 MHz
% derated_dBpcmpMHz = 1.78 / 1.85;  % for liver at 1.8 MHz
% derated_dBpcmpMHz = 1.1 / 1.85;  % for fat at 1.8 MHz
derated_dBpcmpMHz = 0.3;  % standard FDA derating (dB/(cm * MHz))

if(USE_DERATED)
%     derated_dB = derated_dBpcmpMHz .* (freq_ultrasound ./ 1e6) .* (distance .* 100);
    derated_dB = derated_dBpcmpMHz .* 1.8 .* (distance .* 100);
    derated_amplitude_ratio = 10.^(derated_dB./20);
else
    derated_amplitude_ratio = 1;
end
pRms = vRms / (amplifier_amplitude_ratio .* hydrophone_sensitivity .* derated_amplitude_ratio);


intensity_timeAverage = (pRms.^2) ./ (density_water .* speed_sound_water);  % W/m^2

p_amplitude = pRms .* sqrt(2);
MI = (p_amplitude ./ 1e6) ./ sqrt(freq_ultrasound ./ 1e6);  % this assumes that the pressure wave is symmetric about ambient - okay for linear / not high pressure waves


