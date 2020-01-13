%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2020
% All rights reserved.
%========================================


clear
figHandles = findall(0,'Type','figure');
for n = 1:length(figHandles)
    clf(figHandles(n)) % clear figures, but don't close them
end

sweeplabels = {'28_001'};

for sweepidx = 1:length(sweeplabels)
    sweeplabel = sweeplabels{sweepidx};
    sweeplabel = ['BJ_', sweeplabel]
    fprintf('one')
    d = dir(['C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2019-02-16_benchtop_pork_test\matlab\SD_', sweeplabel, '*.mat']);
    load([d.folder, '\', d.name])
    length(dataStore)
    dataStore(1).tp_settings
    for datastoreidx = 1:length(dataStore)
        clf
        time_offset = 0.011529
        for chidx = [3, 4]
            plot(1000*(dataStore(datastoreidx).scope.xData(chidx,:) - time_offset), dataStore(datastoreidx).scope.yData(chidx,:)); hold on
            title(sweeplabel)
        end
        myxlim = xlim;
        xlim([-2, 0.02e3]);
        %ylim([0, 0.09])
        xlabel('time (ms)')
        ylabel('mote voltage (V)')
        title('Mote waveforms')
        drawnow
%         fn_format_and_save_figure(gcf, [pwd filesep 'mote_charging_figures' filesep 'pork_test'], [900.*.60, (550/1).*.75]);
        pause;
    end
    pause;
end




