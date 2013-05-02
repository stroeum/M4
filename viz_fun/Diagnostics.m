clear all
close all
clc

q = 1.60217657e-19; %_C, elementary charge
mi(1) = 5.3135e-26; %_kg, mass of O2+, only single charged ions allowed!
mi(2) = 7.3079e-26; %_kg, mass of CO2+
mi(3) = 2.6567e-26; %_kg, mass of O+
me    = 9.1093e-31; %_kg, mass of e

%% Load data
raw = importdata('../viz_dir/diagnostics.dat');

%% Resize data

% S.N = size(N.data);
t         = raw.data(:,1);
N         = length(t);
Ni(3,N)   = 0.0;
Ne(1,N)   = 0.0;
Fi.s(3,N) = 0.0;
Fi.n(3,N) = 0.0;
Fi.w(3,N) = 0.0;
Fi.e(3,N) = 0.0;
Fi.d(3,N) = 0.0;
Fi.u(3,N) = 0.0;
Fi.s(1,N) = 0.0;
Fi.n(1,N) = 0.0;
Fi.w(1,N) = 0.0;
Fi.e(1,N) = 0.0;
Fi.d(1,N) = 0.0;
Fi.u(1,N) = 0.0;

for l=1:3
    Ni(l,:)   = raw.data(:,2+(l-1)*7);
    Fi.s(l,:) = raw.data(:,3+(l-1)*7);
    Fi.n(l,:) = raw.data(:,4+(l-1)*7);
    Fi.w(l,:) = raw.data(:,5+(l-1)*7);
    Fi.e(l,:) = raw.data(:,6+(l-1)*7);
    Fi.d(l,:) = raw.data(:,7+(l-1)*7);
    Fi.u(l,:) = raw.data(:,8+(l-1)*7);
end
l=4;
Ne(1,:)   = raw.data(:,2+(l-1)*7);
Fe.s(1,:) = raw.data(:,3+(l-1)*7);
Fe.n(1,:) = raw.data(:,4+(l-1)*7);
Fe.w(1,:) = raw.data(:,5+(l-1)*7);
Fe.e(1,:) = raw.data(:,6+(l-1)*7);
Fe.d(1,:) = raw.data(:,7+(l-1)*7);
Fe.u(1,:) = raw.data(:,8+(l-1)*7);

%% Plot
FS=16;
figure(1)
set(gcf,'Units','normalized','OuterPosition',[0 0 .25 1])
subplot(4,1,1)
semilogy(t,Ne,'k-',t,Ni(1,:),'k--',t,Ni(2,:),'k-.',t,Ni(3,:),'k:')
xlabel('t (s)','FontSize',FS);
ylabel('N_\alpha (#)','FontSize',FS);
legend('e','O_2^+','CO_2^+','O^+','Location','best')
legend('Boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(4,1,2)
N = (Ne + sum(Ni,1));
eta = N/N(1);
plot(t,eta,'k-')
ylim([.99 1])
%xlim([0 375])
xlabel('t (s)','FontSize',FS);
ylabel('\Sigma N_\alpha(t) / \Sigma N_\alpha(0)','FontSize',FS);
set(gca,'TickDir','Out','FontSize',FS)

subplot(4,1,3:4)
Fi.S = sum(Fi.s,1);
Fi.N = sum(Fi.n,1);
Fi.W = sum(Fi.w,1);
Fi.E = sum(Fi.e,1);
Fi.D = sum(Fi.d,1);
Fi.U = sum(Fi.u,1);
F.S = Fe.s + Fi.S;
F.N = Fe.n + Fi.N;
F.W = Fe.w + Fi.W;
F.E = Fe.e + Fi.E;
F.D = Fe.d + Fi.D;
F.U = Fe.u + Fi.U;
plot(t,F.S,'k-',t,F.N,'k--',t,F.W,'k-.',t,F.E,'k:',t,F.D,'kv',t,F.U,'k^')
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha(t)','FontSize',FS);
legend('south','north','west','east','down','up','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

figure(2)
set(gcf,'Units','normalized','OuterPosition',[.25 0 .25 1])
subplot(3,2,1)
plot(t,abs(q*Fi.s(1,:)),'r--',t,abs(q*Fi.s(2,:)),'r-.',t,abs(q*Fi.s(3,:)),'r:',t,abs(q*Fe.s(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{South} (A)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,2)
plot(t,abs(q*Fi.n(1,:)),'r--',t,abs(q*Fi.n(2,:)),'r-.',t,abs(q*Fi.n(3,:)),'r:',t,abs(q*Fe.n(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{North} (A)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,3)
plot(t,abs(q*Fi.w(1,:)),'r--',t,abs(q*Fi.w(2,:)),'r-.',t,abs(q*Fi.w(3,:)),'r:',t,abs(q*Fe.w(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{West} (A)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,4)
plot(t,abs(q*Fi.e(1,:)),'r--',t,abs(q*Fi.e(2,:)),'r-.',t,abs(q*Fi.e(3,:)),'r:',t,abs(q*Fe.e(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{East} (A)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,5)
plot(t,abs(q*Fi.d(1,:)),'r--',t,abs(q*Fi.d(2,:)),'r-.',t,abs(q*Fi.d(3,:)),'r:',t,abs(q*Fe.d(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{Down} (A)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,6)
plot(t,abs(q*Fi.u(1,:)),'r--',t,abs(q*Fi.u(2,:)),'r-.',t,abs(q*Fi.u(3,:)),'r:',t,abs(q*Fe.u(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{Up} (A)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

figure(3)
set(gcf,'Units','normalized','OuterPosition',[.5 0 .25 1])
subplot(3,2,1)
plot(t,abs(mi(1)*Fi.s(1,:)),'r--',t,abs(mi(2)*Fi.s(2,:)),'r-.',t,abs(mi(3)*Fi.s(3,:)),'r:',t,abs(me*Fe.s(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{South} (kg/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,2)
plot(t,abs(mi(1)*Fi.n(1,:)),'r--',t,abs(mi(2)*Fi.n(2,:)),'r-.',t,abs(mi(3)*Fi.n(3,:)),'r:',t,abs(me*Fe.n(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{North} (kg/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,3)
plot(t,abs(mi(1)*Fi.w(1,:)),'r--',t,abs(mi(2)*Fi.w(2,:)),'r-.',t,abs(mi(3)*Fi.w(3,:)),'r:',t,abs(me*Fe.w(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{West} (kg/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,4)
plot(t,abs(mi(1)*Fi.e(1,:)),'r--',t,abs(mi(2)*Fi.e(2,:)),'r-.',t,abs(mi(3)*Fi.e(3,:)),'r:',t,abs(me*Fe.e(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{East} (kg/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,5)
plot(t,abs(mi(1)*Fi.d(1,:)),'r--',t,abs(mi(2)*Fi.d(2,:)),'r-.',t,abs(mi(3)*Fi.d(3,:)),'r:',t,abs(me*Fe.d(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{Down} (kg/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,6)
plot(t,abs(mi(1)*Fi.u(1,:)),'r--',t,abs(mi(2)*Fi.u(2,:)),'r-.',t,abs(mi(3)*Fi.u(3,:)),'r:',t,abs(me*Fe.u(1,:)),'b-');
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{Up} (kg/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)


% print(1,'-dps','figure.ps')