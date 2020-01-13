%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2019; Last revision: 2019-02
% All rights reserved.
%========================================




% power harvest at mote calc
% close all
clear
figHandles = findall(0,'Type','figure');
for n = 1:length(figHandles)
    clf(figHandles(n)) % clear figures, but don't close them
end


addpath('C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\scripts');
addpath('C:\Users\david\dkpStore\UCBSF\R\Writing\Resources\figure_colorschemes\DrosteEffect-CubeHelix-00335f8\DrosteEffect-CubeHelix-00335f8')
myColormap = cubehelix(128, 1.72, 1.31, 1.95, 1.30, [0.15, 0.82], [0.15, 0.82]);
numColors = 5;
for n = 1:numColors
    cmapIndex = (floor((n - 1).*((length(myColormap)-1)./(numColors-1)))+1);
    color5{n} = myColormap(cmapIndex, :);
end
set(0,'defaultAxesFontName', 'Arial');
set(0,'defaultTextFontName', 'Arial');
set(groot,'defaultAxesFontName', 'Arial');
set(groot,'defaultTextFontName', 'Arial');





%load('C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2018-02-21_invivo\scope129.mat')

load('C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2019-02-16_benchtop_pork_test\matlab\SD_BJ_35_001__20190216T143731.mat')
figure(301)
%plot(dataStore(1).scope.xData(1,:), dataStore(1).scope.yData(4,:))

v = dataStore(1).scope.yData(4,:)';
stim = dataStore(1).scope.yData(3,:)';
t = dataStore(1).scope.xData(1,:)';
tstart = .00464;
istart = ceil(tstart ./ (t(2)-t(1)));  % 10735
iend = length(v);   %13864

t = t(istart:iend);
t = t - t(1);
%t = t - .09306;
v = v(istart:iend);
v = smooth(v, 360, 'lowess');

figure(1)
plot(t, v, 'LineWidth', 3)
%ylim([0 3.3])

EinCap = .5 .* 4e-6 .* v.^2;

EinCap = smooth(EinCap, 15);
EinCap = smooth(EinCap, 50, 'lowess');

PintoCap = (EinCap(2:end) - EinCap(1:(end-1)))./(t(2:end) - t(1:(end-1)));

% PintoCap = PintoCap(30:(end-30));
PintoCap = smooth(PintoCap, 15);
PintoCap = smooth(PintoCap, 30, 'lowess');

figure(2)
title('EinCap')
plot(t, EinCap)

figure(3)
plot((t(2)-t(1)).*(1:length(PintoCap)), smooth(PintoCap, 30, 'loess'))
title('PintoCap')
ylim([0 3.2e-4])


PintoCapSmoothL = PintoCap(1:509);
PintoCapSmoothR = PintoCap(510:end);

PintoCapSmooth = [smooth(PintoCapSmoothL, 40, 'loess') ; smooth(PintoCapSmoothR, 600, 'loess')];
PintoCapSmooth = smooth(PintoCapSmooth, 40);

figure(4)
plot((t(2)-t(1)).*(1:length(PintoCap)), PintoCapSmooth)
title('PintoCap smooth')
ylim([0 3.2e-4])


figure(5)
plot(v(2:end), PintoCapSmooth)
title('PintoCap smooth')
ylim([0 3.2e-4])



figure(6); % efficiency
Pacoustic = 3.4e-3;  % W
Pacoustic = 2.81e-3;  % W
eta_acoustic_to_Cstore = PintoCapSmooth./Pacoustic;
plot(v(2:end), eta_acoustic_to_Cstore)


figure(1)
xlabel('time (s)')
ylabel('VDD2p5 (V)')
title('mote charting')
set(gca, 'LineWidth', 14)
% fn_format_and_save_figure(1, 'powerup', [900.*.85.*.65, (550/1).*.85.*.60]);




figure(21)
figure(1)
plot(t, v, 'LineWidth', 2)
ylim([0 3.3])
xlabel('Time (s)')
ylabel('VDD (V)')
title('VDD during mote charge-up before stimulation is initiated')
xlim([-.05, 0.35])
ylim([0, 3.3])
drawnow
%fn_format_and_save_figure(gcf, [pwd filesep 'mote_charging_figures' filesep 'mote_charging_v-vs-t'], [900.*.75, (550/1).*.75]);

figure(22)
leftCut = 860;
rightCut = 160;
plot(v((2+leftCut):(end-rightCut)), eta_acoustic_to_Cstore(1+leftCut:(end-rightCut)), 'LineWidth', 2)
% plot(v(2:end), eta_acoustic_to_Cstore)
xlabel('VDD (V)')
ylabel('Efficiency of acoustic to Cstore charge')
title('Efficiency vs V_DD')
xlim([0.5, 3.2])
ylim([0, 0.09])
drawnow
%fn_format_and_save_figure(gcf, [pwd filesep 'mote_charging_figures' filesep 'mote_charging_eta-vs-v'], [900.*.75, (550/1).*.75]);
