clear all
close all
clc

q = 1; %1.60217657e-19; %_C, elementary charge
mi(1) = 5.3135e-26; %_kg, mass of O2+, only single charged ions allowed!
mi(2) = 7.3079e-26; %_kg, mass of CO2+
mi(3) = 2.6567e-26; %_kg, mass of O+
me    = 9.1093e-31; %_kg, mass of e


%% Find folder

cd ..
[upperPath, Folder] = fileparts(pwd);
cd viz_fun


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
Fi.b(3,N) = 0.0;
Fi.t(3,N) = 0.0;
Fe.w(1,N) = 0.0;
Fe.e(1,N) = 0.0;
Fe.s(1,N) = 0.0;
Fe.n(1,N) = 0.0;
Fe.b(1,N) = 0.0;
Fe.t(1,N) = 0.0;

for l=1:3
    Ni(l,:)   = raw.data(:,2+(l-1)*7);
    Fi.w(l,:) = raw.data(:,3+(l-1)*7);
    Fi.e(l,:) = raw.data(:,4+(l-1)*7);
    Fi.s(l,:) = raw.data(:,5+(l-1)*7);
    Fi.n(l,:) = raw.data(:,6+(l-1)*7);
    Fi.b(l,:) = raw.data(:,7+(l-1)*7);
    Fi.t(l,:) = raw.data(:,8+(l-1)*7);
end
l=4;
Ne(1,:)   = raw.data(:,2+(l-1)*7);
Fe.w(1,:) = raw.data(:,3+(l-1)*7);
Fe.e(1,:) = raw.data(:,4+(l-1)*7);
Fe.s(1,:) = raw.data(:,5+(l-1)*7);
Fe.n(1,:) = raw.data(:,6+(l-1)*7);
Fe.b(1,:) = raw.data(:,7+(l-1)*7);
Fe.t(1,:) = raw.data(:,8+(l-1)*7);

%% Plot
FS=16;
figure(1)
set(gcf,'Units','normalized','OuterPosition',[0 0 .25 1],'Color',[1 1 1])

subplot(4,1,1)
semilogy(t,Ne,'b-',t,Ni(1,:),'r--',t,Ni(2,:),'r-.',t,Ni(3,:),'r:')
xlabel('t (s)','FontSize',FS);
ylabel('N_\alpha (#)','FontSize',FS);
legend('e','O_2^+','CO_2^+','O^+','Location','best')
legend('Boxoff')
axis tight
set(gca,'TickDir','Out','FontSize',FS)
title(Folder)

subplot(4,1,2)
N = (Ne + sum(Ni,1));
eta = N/N(1);
plot(t,eta,'k-')
axis tight
%ylim([.99 1])
%xlim([0 375])
xlabel('t (s)','FontSize',FS);
ylabel('\Sigma_\alpha N_\alpha(t) / \Sigma N_\alpha(0)','FontSize',FS);
set(gca,'TickDir','Out','FontSize',FS)

subplot(4,1,3:4)
% Sum_i ni vi %
Fi.S = sum(Fi.s,1);
Fi.N = sum(Fi.n,1);
Fi.W = sum(Fi.w,1);
Fi.E = sum(Fi.e,1);
Fi.B = sum(Fi.b,1);
Fi.T = sum(Fi.t,1);
% en masse flux  Sum_alpha n_alpha v_alpha = ne ve + Sum_i ni vi %
F.S = Fe.s + Fi.S;
F.N = Fe.n + Fi.N;
F.W = Fe.w + Fi.W;
F.E = Fe.e + Fi.E;
F.B = Fe.b + Fi.B;
F.T = Fe.t + Fi.T;
% charged flux  Sum_alpha n_alpha v_alpha = -ne ve + Sum_i ni vi %
% F.S = -Fe.s + Fi.S;
% F.N = -Fe.n + Fi.N;
% F.W = -Fe.w + Fi.W;
% F.E = -Fe.e + Fi.E;
% F.B = -Fe.b + Fi.B;
% F.T = -Fe.t + Fi.T;

plot(t,F.S,'k-',t,F.N,'k--',t,F.W,'k-.',t,F.E,'k:',t,F.B,'kv',t,F.T,'k^')
%plot(t,F.S+F.N+F.W+F.E+F.B+F.T,'k')
xlabel('t (s)','FontSize',FS);
ylabel('F^\beta=(\Sigma_in_iv_i^\beta+n_ev_e^\beta).dS (#/s)','FontSize',FS);
axis tight
legend('south','north','west','east','down','up','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

figure(2)
set(gcf,'Units','normalized','OuterPosition',[.25 0 .25 1],'Color',[1 1 1])

minF = -15e22;
maxF = +15e22;

subplot(3,2,1)
plot(t,Fi.s(1,:),'r--',t,Fi.s(2,:),'r-.',t,Fi.s(3,:),'r:',t,Fe.s(1,:),'b-');
ylim([minF maxF])
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{South} (#/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)
title(Folder)

subplot(3,2,2)
plot(t,Fi.n(1,:),'r--',t,Fi.n(2,:),'r-.',t,Fi.n(3,:),'r:',t,Fe.n(1,:),'b-');
ylim([minF maxF])
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{North} (#/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,3)
plot(t,Fi.w(1,:),'r--',t,Fi.w(2,:),'r-.',t,Fi.w(3,:),'r:',t,Fe.w(1,:),'b-');
ylim([minF maxF])
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{West} (#/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,4)
plot(t,Fi.e(1,:),'r--',t,Fi.e(2,:),'r-.',t,Fi.e(3,:),'r:',t,Fe.e(1,:),'b-');
ylim([minF maxF])
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{East} (#/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,5)
plot(t,Fi.b(1,:),'r--',t,Fi.b(2,:),'r-.',t,Fi.b(3,:),'r:',t,Fe.b(1,:),'b-');
ylim([minF maxF])
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{Down} (#/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

subplot(3,2,6)
plot(t,Fi.t(1,:),'r--',t,Fi.t(2,:),'r-.',t,Fi.t(3,:),'r:',t,Fe.t(1,:),'b-');
ylim([minF maxF])
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{Up} (#/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)

% figure(3)
% set(gcf,'Units','normalized','OuterPosition',[.5 0 .25 1])
% subplot(3,2,1)
% plot(t,mi(1)*Fi.s(1,:),'r--',t,mi(2)*Fi.s(2,:),'r-.',t,mi(3)*Fi.s(3,:),'r:',t,me*Fe.s(1,:),'b-');
% xlabel('t (s)','FontSize',FS);
% ylabel('F_\alpha^{South} (kg/s)','FontSize',FS);
% legend('O_2^+','CO_2^+','O^+','e','location','best')
% legend('boxoff')
% set(gca,'TickDir','Out','FontSize',FS)
% 
% subplot(3,2,2)
% plot(t,mi(1)*Fi.n(1,:),'r--',t,mi(2)*Fi.n(2,:),'r-.',t,mi(3)*Fi.n(3,:),'r:',t,me*Fe.n(1,:),'b-');
% xlabel('t (s)','FontSize',FS);
% ylabel('F_\alpha^{North} (kg/s)','FontSize',FS);
% legend('O_2^+','CO_2^+','O^+','e','location','best')
% legend('boxoff')
% set(gca,'TickDir','Out','FontSize',FS)
% 
% subplot(3,2,3)
% plot(t,mi(1)*Fi.w(1,:),'r--',t,mi(2)*Fi.w(2,:),'r-.',t,mi(3)*Fi.w(3,:),'r:',t,me*Fe.w(1,:),'b-');
% xlabel('t (s)','FontSize',FS);
% ylabel('F_\alpha^{West} (kg/s)','FontSize',FS);
% legend('O_2^+','CO_2^+','O^+','e','location','best')
% legend('boxoff')
% set(gca,'TickDir','Out','FontSize',FS)
% 
% subplot(3,2,4)
% plot(t,mi(1)*Fi.e(1,:),'r--',t,mi(2)*Fi.e(2,:),'r-.',t,mi(3)*Fi.e(3,:),'r:',t,me*Fe.e(1,:),'b-');
% xlabel('t (s)','FontSize',FS);
% ylabel('F_\alpha^{East} (kg/s)','FontSize',FS);
% legend('O_2^+','CO_2^+','O^+','e','location','best')
% legend('boxoff')
% set(gca,'TickDir','Out','FontSize',FS)
% 
% subplot(3,2,5)
% plot(t,mi(1)*Fi.b(1,:),'r--',t,mi(2)*Fi.b(2,:),'r-.',t,mi(3)*Fi.b(3,:),'r:',t,me*Fe.b(1,:),'b-');
% xlabel('t (s)','FontSize',FS);
% ylabel('F_\alpha^{Down} (kg/s)','FontSize',FS);
% legend('O_2^+','CO_2^+','O^+','e','location','best')
% legend('boxoff')
% set(gca,'TickDir','Out','FontSize',FS)
% 
% subplot(3,2,6)
% plot(t,mi(1)*Fi.t(1,:),'r--',t,mi(2)*Fi.t(2,:),'r-.',t,mi(3)*Fi.t(3,:),'r:',t,me*Fe.t(1,:),'b-');
% xlabel('t (s)','FontSize',FS);
% ylabel('F_\alpha^{Up} (kg/s)','FontSize',FS);
% legend('O_2^+','CO_2^+','O^+','e','location','best')
% legend('boxoff')
% set(gca,'TickDir','Out','FontSize',FS)


print(1,'-dps',[Folder,'_TotalFlowPerBoundary.ps'])
print(2,'-dps',[Folder,'_FlowPerSpeciesPerBoundary.ps'])