%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2018-06-20
% All rights reserved.
%========================================


close all

open('C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\scripts\figures_out_bmod_norm\C98_small25\xtransHydrophoneLarge25Intensity.fig')  % small 25 (derated?)
% open('C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\scripts\figures_out_bmod_norm\C68_large25\xtransHydrophoneLarge25Intensity.fig')  % large 25 (not derated?)

h = gcf;

axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');

% FOR LARGE, CHANGE TO NOT-DERATED FROM DERATED AND CHANGE TO 28.9 V FROM 25.6 V
ydata = ydata .* (1006 ./ 547.3) .* (28.9 ./ 25.6).^2;

% xdata is in microns; y data is in mW/cm^2
x = xdata ./ 1e6;  % now xdata is in m
y = ydata .* 10000 ./ 1000;  % W/m^2

figure(101)
plot(x,y)
ylim([-200 16000])
figure(102)


% average accross y-axis to get radial distribution
ptsToUse = 98;  % 98  40
zeroIdx = find(x == 0);
xavg = x((zeroIdx):-1:(zeroIdx-ptsToUse));
yavg = [y(zeroIdx) mean([y((zeroIdx+1):1:(zeroIdx+ptsToUse)); y((zeroIdx-1):-1:(zeroIdx-ptsToUse))], 1)];
plot(xavg,yavg)

for ring = 1:(length(xavg)-1)
    A(ring) = pi * (xavg(ring+1).^2 - xavg(ring).^2);
    I(ring) = mean([yavg(ring+1), yavg(ring)]);
    P(ring) = A(ring).*I(ring);
    Pencircled(ring) = sum(P(1:ring));
end

figure(103)
plot(xavg(2:end),Pencircled)




areaSquare = 750e-6 .^ 2;
rCircle = (areaSquare ./ pi).^(.5);

vq = interp1(xavg(2:end), Pencircled, rCircle)

