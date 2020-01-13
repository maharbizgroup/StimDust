%========================================
% StimDust
% Author: David K. Piech
% University of California, Berkeley
% email address: piech@berkeley.edu
% Website: 
%     https://maharbizgroup.wordpress.com/
%     http://carmenalab.org/
%     https://people.eecs.berkeley.edu/~rikky/Home.html
% 2018; Last revision: 2018-06-16
% All rights reserved.
%========================================



clear
figHandles = findall(0,'Type','figure');
for n = 1:length(figHandles)
    clf(figHandles(n)) % clear figures, but don't close them
end

%================================================================
%benchtop
load('D:\yesfh\stimDustData\backscatter_and_scope_data\scope_D3.mat')
figure(10); hold on
plot(xData(2,:), yData(2,:))
plot(xData(1,:), yData(1,:))
plot(xData(3,:), yData(3,:))
plot(xData(4,:), yData(4,:))

t = xData(4,:);
t = t - 142.365e-6 + 84.78e-8;
v = yData(4,:);
v = v;% + .025;
v = smooth(v, 90);
idxStart = find(abs(t) < 5e-9);
idxStart = idxStart(end)
idxEnd = find(abs(t - 7.364e-5)<1e-8);
idxEnd = idxEnd(end)-7 + 600;

t = t(idxStart:idxEnd);
v = v(idxStart:idxEnd);


R = 3000;  % Ohm
C = 22e-9;  % Farad

%========================================================================
%invivo
M = csvread('C:\Users\david\dkpStore\UCBSF\R\ND\StimDust\experiment\2018-02-21_invivo\scope\scope_127.csv', 2, 0);
t = M(:,1);
v = M(:,4);
v = smooth(v, 90);

idxstart = 26000;
idxend = 42000;
t = t(idxstart:idxend);
v = v(idxstart:idxend);

R = 4400;  % Ohm
C = 17.9e-9;  % Farad

%========================================================================






figure(1); clf; hold on
plot(t,v)
% ylim([0 2])


Dt = t(3) - t(2);

% vC = Q/C = i Dt / C
% vR = IR

vdd = 1.86;


iStart = 386e-6;  % A
iEnd = 193e-6;  % A


ihat = iStart .* ones(length(t), 1);

idend = length(t);
idcomp = 3408;
ihat(idcomp:end) = linspace(iStart,iEnd, (length(t) - idcomp + 1));  % iStart:(-(iEnd-iStart)./(idend-idcomp)):iEnd;

vRhat = R .* ihat;
DvChat = ihat .* Dt ./ C;
vChat = cumsum(DvChat);
vhat = vRhat + vChat;


% plot(t, vRhat)
% plot(t, vChat)
% plot(t, vhat)


% figure(2); clf; hold on
% PR = vRhat.^2 ./ R;
% ER = cumsum(PR).*Dt;
% EC = .5 .* C .* vChat.^2;
% plot(t, ER)
% plot(t, EC)


for n = 1:length(t)
    if n == 1
        ih(n) = (v(n) - 0) ./ (R + Dt./C);
        vCh(n) = 0;
    else
        ih(n) = (v(n) - vCh(n-1)) ./ (R + Dt./C);
        vCh(n) = vCh(n-1) + ih(n) .* Dt ./ C;
    end
    
    vRh(n) = ih(n) .* R;
    
    PR(n) = ih(n).^2 .* R;
    EtotR(n) = sum(PR(1:n)) .* Dt;
    EtotC(n) = .5 .* C .* vCh(n) .^ 2;
    Etot(n) = EtotR(n) + EtotC(n);
    
    P_iv(n) = ih(n) .* v(n);
    Etot_iv(n) = sum(P_iv(1:n)) .* Dt;
    
end
Q = cumsum(ih) .* Dt;

figure(3); hold on
% plot(t, ih)
plot(t, vCh)
plot(t, vRh)

figure(4); hold on
plot(t, ih)
c
figure(5); hold on
plot(t, Q)

figure(6); hold on
% plot(t, EtotR)
% plot(t, EtotC)
plot(t, Etot_iv)
plot(t, Etot, 'b.')