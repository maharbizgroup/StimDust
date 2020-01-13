function [value_out] = fn_txdr_voltage_scaling(txdr_name, p_or_i, value)
%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2019-05-13
% All rights reserved.
%========================================



if strcmp(txdr_name, 'large')
    vTx = [5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 32];
    p_peakRare = 1.0e+09 * [0.1930, 0.3408, 0.4854, 0.6359, 0.7823, 0.9097, 1.0322, 1.1513, 1.2681, 1.3554, 1.4354, 1.5048, 1.5647, 1.6019, 1.6134];
    intensity_timeAverage = 1.0e+11 * [0.1209, 0.3771, 0.7650, 1.3127, 1.9869, 2.6871, 3.4595, 4.3035, 5.2208, 5.9649, 6.6897, 7.3524, 7.9487, 8.3316, 8.4512];
else
    vTx = [5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 32];
    p_peakRare = 1.0e+09 * [0.2309, 0.3699, 0.4998, 0.6332, 0.7631, 0.8894, 1.0152, 1.1344, 1.2532, 1.3687, 1.4765, 1.5916, 1.6907, 1.7886, 1.8246];
    intensity_timeAverage = 1.0e+12 * [0.0173, 0.0444, 0.0811, 0.1302, 0.1891, 0.2568, 0.3346, 0.4178, 0.5099, 0.6082, 0.7078, 0.8224, 0.9280, 1.0386, 1.0809];
end

if strcmp(p_or_i, 'p') | strcmp(p_or_i, 'v')
    y = p_peakRare;
else
    y = intensity_timeAverage;
end

value_out = interp1(vTx, y, value);