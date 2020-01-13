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

figHandles = findall(0,'Type','figure'); for n = 1:length(figHandles); clf(figHandles(n)); end   % clear figures, but don't close them

if ~ exist('data', 'var')
    load('C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\scripts\emg_process_dataout\pedot1\dataWorking5.mat');
end


% set colormap
% note: need cubehelix function from mathworks
addpath('C:\Users\david\dkpStore\UCBSF\R\Writing\Resources\figure_colorschemes\DrosteEffect-CubeHelix-00335f8\DrosteEffect-CubeHelix-00335f8')
myColormap = cubehelix(1024, 1.72, 1.31, 1.95, 1.30, [0.15, 0.82], [0.15, 0.82]);
numColors = 5;
for n = 1:numColors
    cmapIndex = (floor((n - 1).*((length(myColormap)-1)./(numColors-1)))+1);
    color5{n} = myColormap(cmapIndex, :);
    colors = color5;
end

myColormap = myColormap(floor(0.13*length(myColormap)):floor(0.95*length(myColormap)), :);

set(0,'defaultAxesFontName', 'Arial');
set(0,'defaultTextFontName', 'Arial');
set(groot,'defaultAxesFontName', 'Arial');
set(groot,'defaultTextFontName', 'Arial');


figure(1); hold on
device_num = 2;
ampval = 7;

nominal_stim_durations = [400, 300, 200, 100, 400, 300, 200, 100];
is_plateds = [0, 0, 0, 0, 1, 1, 1, 1];
HforLegend = [];

for n = 1:length(nominal_stim_durations)
    [stims, line] = grab_from_database(data, ampval, nominal_stim_durations(n), device_num, is_plateds(n));
    size(stims, 1)
    
    
    % deal with 
%     if nominal_stim_durations(n) == 200 && is_plateds(n) == 0
%         figure(10+n)
%         for k = 1:size(stims, 1)
%             plot(stims(k,:)); hold on
%             pause
%         end
%     end
%     figure(10+n)
%     plot(stims')
%     figure(1)

    %colorThisCondition = myColormap(floor((1 ./ 4) .* length(myColormap)), :);
    colorThisCondition = myColormap(floor(((((nominal_stim_durations(n) ./ 100) - 1) * 0.999 ./ 3) .* length(myColormap)) + 1), :);
    c_hsv = rgb2hsv(colorThisCondition);
    if is_plateds(n)
        c_hsv(2) = 1.3.*c_hsv(2);
        if c_hsv(2) > 1
            c_hsv(2) = 1;
        end
    else
        c_hsv(2) = 0.4.*c_hsv(2);
    end
    colorThisCondition = hsv2rgb(c_hsv);

    H = shadedErrorBar(line.vstimTime .* 1e6, mean(stims, 1), std(stims, 0, 1), 'lineProps', {'-k', 'Color', colorThisCondition, 'LineWidth', 1.00, 'DisplayName', sprintf('%3.0f %cs N=%u', line.duration .* 1e6, 956, size(stims, 1))});
    
    HforLegend(n) = H.mainLine;
    ax = get(gca); ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
    title([sprintf('Electrodes voltage'), sprintf('pedot_device%u_ampval%u', device_num, ampval)], 'Interpreter', 'none')
    %xlim(emgWindowView .* 1000)
    xlabel('time (\mus)'); ylabel('V_{electrodes} (V)');
    yLimCurr = ylim;
    ylim([-max(abs(yLimCurr)), max(abs(yLimCurr))])
    ylim([-.25, 2.5])

    figure(1)
    lgd1 = legend('show', 'Location', 'northeast');
    lgd1 = legend(HforLegend);
    legend boxoff;
    lgd1.FontSize = 12; lgd1.EdgeColor = [.65, .65, .65];
    
    % export to csv
    if nominal_stim_durations(n) == 400
%         disp('here')
        vmean = mean(stims, 1)';
        tstart = find(vmean > 0.05);
        tstart = tstart(1);
        t = line.vstimTime' - line.vstimTime(tstart);
%         csvwrite(sprintf('device%u_ampval%u.csv', device_num, ampval), [t, vmean]);
%         pause
        
        %csvwrite(sprintf('animalD1_ampval0.csv', device_num, ampval), [1e-6.*dataObjs(1).XData', dataObjs(1).YData']);
        %csvwrite('device_1pre.txt', [f', (m_3pre' .* sind(p_3pre')), (m_3pre' .* cosd(p_3pre'))]);
        %csvwrite('device_1pre.txt', [f', (1e6*m_3pre' .* cosd(p_3pre')), (1e6.*m_3pre' .* sind(p_3pre'))]);
    end
        
    
end

% fn_format_and_save_figure(1, [pwd filesep 'emg_process_dataout' filesep 'pedot1' filesep sprintf('pedot_device%u_ampval%u', device_num, ampval)], [900.*.65 550.*.65])  % [width height]

