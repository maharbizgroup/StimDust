%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2020; Last revision: 2020-01-13
% All rights reserved.
%========================================



%=========================================================================
% initialize
clear
clc
figHandles = findall(0,'Type','figure');
for n = 1:length(figHandles)
    clf(figHandles(n)) % clear figures, but don't close them
end



%=========================================================================
% metadata

data = struct('lineNum', {});
data(1).ampval = 10;
data(2).ampval = 7;
data(3).ampval = 6;
data(4).ampval = 5;
data(5).ampval = 4;
data(6).ampval = 3;
data(7).ampval = 2;
data(8).ampval = 1;
data(9).ampval = 0;

data(1).scopenum = '22';
data(2).scopenum = '23';
data(3).scopenum = '24';
data(4).scopenum = '26';
data(5).scopenum = '25';
data(6).scopenum = '27';
data(7).scopenum = '28';
data(8).scopenum = '29';
data(9).scopenum = '30';


data(1).current = 400e-6;
data(2).current = 400e-6;
data(3).current = 350e-6;
data(4).current = 250e-6;
data(5).current = 300e-6;
data(6).current = 150e-6;
data(7).current = 100e-6;
data(8).current = 50e-6;
data(9).current = 50e-6;

% use every fourth slice because EMG capture was every 0.25 seconds, but stimulator PRF was 1 Hz.
blank = 3;
interval = 4;
count = 55;
data(2).slices_to_use = (blank+1):interval:((count-1)*interval+blank); % 3:4:55*4+3
blank = 0;
interval = 4;
count = 25;
data(3).slices_to_use = (blank+1):interval:((count-1)*interval+blank);
blank = 0;
interval = 4;
count = 42;
data(4).slices_to_use = (blank+1):interval:((count-1)*interval+blank);
blank = 1;
interval = 4;
count = 52;
data(5).slices_to_use = (blank+1):interval:((count-1)*interval+blank);
blank = 2;
interval = 4;
count = 30;
data(6).slices_to_use = (blank+1):interval:((count-1)*interval+blank);
blank = 2;
interval = 4;
count = 28;
data(7).slices_to_use = (blank+1):interval:((count-1)*interval+blank);
blank = 0;
interval = 4;
count = 43;
data(8).slices_to_use = (blank+1):interval:((count-1)*interval+blank);
blank = 1;
interval = 4;
count = 37;
data(9).slices_to_use = (blank+1):interval:((count-1)*interval+blank);



lineNums_to_process = 2:9;

for lineNum = lineNums_to_process
    data(lineNum).breakout_num = num2str(lineNum);
end

%=========================================================================
% v stim

figure(1); hold on
for lineNum = lineNums_to_process
    data(lineNum).v_electrodes = ...
        csvread(['C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2017-07-27_in-vivo_breakout\Vstim_scope_data\', ...
        'scope_', data(lineNum).scopenum, '.csv'], 22, 0);
    plot(data(lineNum).v_electrodes(:,1), data(lineNum).v_electrodes(:,2))
end

% create vstim and vdd low-pass filter and apply filter
firNumTaps = 900;
cutoffFreq = 1000000;   % Hz
scopeSamplingFreq = 1./(data(lineNum).v_electrodes(2,1) - data(lineNum).v_electrodes(1,1));
nyquistRate = scopeSamplingFreq./2;  
Wn = cutoffFreq./nyquistRate;
firCoefsScopeLowPass1 = fn_gaussfiltcoef(scopeSamplingFreq, cutoffFreq);

figure(2); hold on
for lineNum = lineNums_to_process
    data(lineNum).v_electrodes_filtered = filtfilt(firCoefsScopeLowPass1, 1, data(lineNum).v_electrodes(:,2));
    plot(data(lineNum).v_electrodes(:,1), data(lineNum).v_electrodes_filtered, 'DisplayName', data(lineNum).breakout_num)
    legend
end


%=========================================================================
% emg

emg_gain = 100;

for lineNum = lineNums_to_process
%     emg_load = load(['C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2017-07-27_in-vivo_breakout\breakout_', data(lineNum).breakout_num, '.txt']);
%     save(['C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2017-07-27_in-vivo_breakout\breakout_', data(lineNum).breakout_num, '.mat'], 'emg_load');
%     continue
    load(['C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2017-07-27_in-vivo_breakout\breakout_', data(lineNum).breakout_num, '.mat']);
    data(lineNum).emgdata_all = emg_load;
    data(lineNum).emgdata_all(:,2) = data(lineNum).emgdata_all(:,2) ./ emg_gain;
    indices_zero_time = find(data(lineNum).emgdata_all(:,1) == 0);
    num_sample_per_emg_capture = indices_zero_time(2) - indices_zero_time(1) - 1;
    data(lineNum).emg_sample_rate = 1 ./ (data(lineNum).emgdata_all(2,1) - data(lineNum).emgdata_all(1,1));
    for n = 1:length(indices_zero_time)-1
        data(lineNum).emgdata_slice{n} = data(lineNum).emgdata_all(indices_zero_time(n):(indices_zero_time(n)+num_sample_per_emg_capture),2);
    end
    emgdata_slice_select = cell(length(data(lineNum).slices_to_use));
    for ii = 1:length(data(lineNum).slices_to_use)
        emgdata_slice_select{ii} = data(lineNum).emgdata_slice{data(lineNum).slices_to_use(ii)};
    end
    data(lineNum).emgdata_slice = emgdata_slice_select;

%     figure(20+lineNum); hold on
    figure(20); hold on
    subplot(1, length(lineNums_to_process), find(lineNums_to_process == lineNum)); hold on
    title(num2str(data(lineNum).current))
    for n = 1:floor(1 * length(data(lineNum).emgdata_slice))
        plot(data(lineNum).emgdata_slice{n})
    end
    
    for n = 1:length(data(lineNum).emgdata_slice)
        data(lineNum).emg_amp_cell{n} = abs(min(data(lineNum).emgdata_slice{n}) - mean(data(lineNum).emgdata_slice{n}));
    end
    data(lineNum).emg_amp = cell2mat(data(lineNum).emg_amp_cell);
    
%     plot(data(lineNum).emgdata_all(:,1), data(lineNum).emgdata_all(:,2))
end


emgSamplingFreq = data(lineNum).emg_sample_rate;
firNumTaps = 698;
cutoffFreq = 50;   % Hz
nyquistRate = emgSamplingFreq./2;  
Wn = [cutoffFreq./nyquistRate];
firCoefsRawHighPassFindTrigger = fir1(firNumTaps, Wn, 'high');

firNumTaps = 80;
cutoffFreq = 2550;   % Hz
nyquistRate = emgSamplingFreq./2;  
Wn = [cutoffFreq./nyquistRate];
firCoefsRawLowPassFindTrigger = fir1(firNumTaps, Wn, 'low');


for lineNum = lineNums_to_process
    for n = 1:length(data(lineNum).emgdata_slice)
        data(lineNum).emgdata_slice_fitered{n} = filtfilt(firCoefsRawHighPassFindTrigger, 1, data(lineNum).emgdata_slice{n});
    end

%     figure(40+lineNum); hold on
    figure(40); hold on
    subplot(1, length(lineNums_to_process), find(lineNums_to_process == lineNum)); hold on
    title(num2str(data(lineNum).current))
    for n = 1:floor(1 * length(data(lineNum).emgdata_slice_fitered))
        plot(data(lineNum).emgdata_slice_fitered{n})
    end
    
    
%     figure(50+lineNum); hold on
%     for n = 1:floor(1 * length(data(lineNum).emgdata_slice_fitered))
%         subplot(4, ceil(floor(1 * length(data(lineNum).emgdata_slice_fitered))/4), n); hold on
%         plot(data(lineNum).emgdata_slice_fitered{n})
%         ylim([-5e-4, 5e-4])
%     end

    

    % emg amplitude single-sided (just above baseline)
    data(lineNum).emgAmplitudesSingle = zeros(length(data(lineNum).emgdata_slice_fitered), 1);
    for n = 1:floor(1 * length(data(lineNum).emgdata_slice_fitered))
        data(lineNum).emgAmplitudesSingle(n) = abs(min(data(lineNum).emgdata_slice_fitered{n}));
    end
    data(lineNum).emgAmplitudeSingleMean = mean(data(lineNum).emgAmplitudesSingle);
    data(lineNum).emgAmplitudeSingleStd = std(data(lineNum).emgAmplitudesSingle);

end

% 
% currents = [];
% amps = [];
% for lineNum = lineNums_to_process
%     currents(end+1) = data(lineNum).current;
%     amps(end+1) = mean(data(lineNum).emg_amp);
% end
% 
% figure(101); hold on
% plot(currents, amps)






figure(4); clf; figure(5); clf; figure(6); clf; figure(7); clf; figure(8); clf; figure(9); clf;
recruitAmpPopX = [];
recruitAmpPopY = [];

%======== time-domain plot and collect data
for n = 1:length(lineNums_to_process)
    lineNum = lineNums_to_process(n);

    %======== X data
    recruitAmpX(n) = data(lineNum).current;

    %======== Y and error data
    recruitAmpY(n) = data(lineNum).emgAmplitudeSingleMean;
    recruitAmpErr(n) = data(lineNum).emgAmplitudeSingleStd;

    %======== population data
    recruitAmpPopY = [recruitAmpPopY; data(lineNum).emgAmplitudesSingle];
    recruitAmpPopX = [recruitAmpPopX; (ones((length(recruitAmpPopY) - length(recruitAmpPopX)), 1).*data(lineNum).current)];
end





%======== recruitment curve sigmoid fit
[xData, yData] = prepareCurveData( recruitAmpPopX, recruitAmpPopY );

% Set up fittype and options.
ft = fittype( 'a/(1+exp(-b*(x-c))) + d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMinChange = 1e-12;
opts.Display = 'Off';
opts.Lower = [0 -Inf -Inf 0];
opts.MaxFunEvals = 60000;
opts.MaxIter = 40000;
opts.StartPoint = [0.004 50000 0.0001 1e-06];
opts.TolFun = 1e-12;
opts.TolX = 1e-08;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

x = -1e-3:.1e-6:1e-3;
yFit = fitresult.a./(1+exp(-fitresult.b.*(x - fitresult.c))) + fitresult.d;


%======== plot recruitment curve
colorLineSweeps = [0, 0.1563, 1.0000];

figure(8); clf; hold on;
plot(x .* 1e6, yFit .* 1e3, 'LineWidth', 1.25);
errorbar(recruitAmpX .* 1e6, recruitAmpY .* 1e3, recruitAmpErr .* 1e3, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', .1, 'CapSize', 10, 'LineWidth', 1.25)  % data in s and V; plot in us and mV
plot(recruitAmpPopX .* 1e6, recruitAmpPopY .* 1e3, 'b.', 'MarkerSize', 6)
ax = get(gca); ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
ylabel('emg amplitude (mV)')
xlim([(min(recruitAmpX .* 1e6) - 50) (max(recruitAmpX .* 1e6) + 50)])
ylim(1000.*[-.0006, .006])

figure(9); hold off
H9 = shadedErrorBar(recruitAmpX .* 1e6, recruitAmpY .* 1e3, recruitAmpErr .* 1e3, 'lineProps', {'-k', 'Color', colorLineSweeps, 'DisplayName', sprintf('%3.0f %cA', data(lineNum).current .* 1e6, 956)});
ax = get(gca); ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
ylabel('emg amplitude (mV)')
ylim(1000.*[-.001, .005])

figure(8); xlabel('stim current (\muA)'); title('emg amplitude vs stim current'); xlim(1e6 .* [0 450e-6])
figure(9); xlabel('stim current (\muA)'); title('emg amplitude vs stim current'); xlim(1e6 .* [0 450e-6])


%======== save figures
xlabel('X Axis', 'FontSize',12)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)

figuresFoldername = 'C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\scripts\figures_out_pilot_data_invivo_breakout_process_2019-12-30';
fn_format_and_save_figure(8, [pwd filesep figuresFoldername filesep 'emgamp_sweep'], [0])

clear H;
close 7; close 9;  % close figures we don't care about


