function [out_matrix, out_line] = grab_from_database(data, ampval, nominal_stim_duration, device_num, is_plated)
%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2019
% All rights reserved.
%========================================

out_matrix = [];
out_line = [];
for lineNum = 1:length(data)
    if (data(lineNum).ampval == ampval) && (data(lineNum).nominal_stim_duration == nominal_stim_duration) && (data(lineNum).device_num == device_num) && (data(lineNum).is_plated == is_plated)
        out_matrix = [out_matrix; data(lineNum).vstim];
        if isempty(out_line)
            out_line = data(lineNum);
            
        end
    end
end