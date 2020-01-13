%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2019-05
% All rights reserved.
%========================================


close all
figHandles = findall(0,'Type','figure');
for n = 1:length(figHandles)
    clf(figHandles(n)) % clear figures, but don't close them
end

figures_folder = 'C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\scripts\2019-05-12_vpiezo_angle\';


uiopen('C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\scripts\Vpiezo.fig',1);
% fig = gcf;
% dataObjs = findobj(fig, '-property', XData');
% x1 = dataObjs(1).YData;
% dataObjs = findobj(fig, '-property', YData');
% y1 = dataObjs(1).YData;
% dataObjs = findobj(fig, '-property', ZData');
% z1 = dataObjs(1).YData;


ch = get(gca, 'ch');
X = get(ch, 'XData');  % translation
Angle = get(ch, 'YData');  % angle
Z = get(ch, 'ZData');  % Vpkpk piezo

x = X(1,:);
angle = Angle(:,1);

x0_idx = 19;
angle0_idx = 9;

x0 = x;
angle0 = angle;
z0 = Z(:, x0_idx);


figure; hold on
x_idx_offset = 6;

zd = [];

min_power_at_which_mote_can_operate = 56;

count = 3;
for x_idx = (-x_idx_offset:1:x_idx_offset) + x0_idx
    z = Z(:, x_idx);
    z_out = z./max(z);
    z_out = min_power_at_which_mote_can_operate./z_out.^2;
    zd(count,:) = z_out;
    plot(angle, z_out, 'LineWidth', 1.5, 'Color', [.25, .25, .25, .15])
    count = count+1;
end
xlim([-65, 65])
% ylim([0, 1.1])
ylim([0, 800])
xlim_here = xlim()
plot(xlim_here, [720, 720], 'g')
title('normalized acoustic mote piezo harvest vpkpk versus angle')
xlabel('angle offset from optimal')
ylabel('normalized pressure amplitude)')
% fn_format_and_save_figure(gcf, [figures_folder 'angle_sensitivity_acoustic_empirical' num2str(x_idx_offset)], [900*.75 550*.75])
% fn_format_and_save_figure(gcf, [figures_folder 'angle_sensitivity_acoustic_power_necessary_56_and_FDA' num2str(x_idx_offset)], [900*.75 550*.75])
set(gca,'XMinorTick','on','YMinorTick','on')
ylim([0, 200])
[x_click,y_click] = ginput(2)

figure; hold on
plot(x, (Z(angle0_idx, :)).^2)
plot([x(-x_idx_offset+x0_idx), x(-x_idx_offset+x0_idx)], ylim(), 'Color', [.7, .7, .7])
plot([x(+x_idx_offset+x0_idx), x(+x_idx_offset+x0_idx)], ylim(), 'Color', [.7, .7, .7])


angle = angle-1;

zds = sort(zd, 1);
z_use = zds(ceil(0.25*size(zds,1)),:);
figure; hold on
plot(angle, mean(zd,1))
plot(angle, z_use)


angle_interp = -min(abs([angle(1), angle(end)])):1:min(abs([angle(1), angle(end)]));
angle_interp = -25:1:25;
z_interp = interp1(angle, z_use, angle_interp);

plot(angle_interp, z_interp)

angle_fold = angle_interp(ceil(length(angle_interp)./2):length(angle_interp));
z_fold = [z_interp(ceil(length(z_interp)./2):length(z_interp));
          z_interp(ceil(length(z_interp)./2):-1:1)];

plot(angle_fold, z_fold(1,:))
plot(angle_fold, z_fold(2,:))

angle_out = angle_fold;
v_out = min(z_fold, [], 1);
p_out = v_out.^2;

figure; hold on
plot(angle_out, v_out)
plot(angle_out, p_out)


a = angle_out;
v = v_out;
p = p_out;

% save('p_v_vs_angle.mat', 'a', 'v', 'p')


