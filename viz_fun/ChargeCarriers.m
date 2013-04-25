clear all
close all
clc
%% Load data
N.data = load('../viz_dir/CarriersCounts.dat');
t = load('../output/t.dat');

%% Resize data

S.N = size(N.data);
N.e    = N.data(:,1);
N.O2p  = N.data(:,2);
N.CO2p = N.data(:,3);
N.Op   = N.data(:,4);
S.t = length(t);
if (S.N(1)<S.t)
    t = t(1:S.N(1));
else
    N.e    = N.e(1:S.t);
    N.O2p  = N.O2p(1:S.t);
    N.CO2p = N.CO2p(1:S.t);
    N.Op   = N.Op(1:S.t);
end

%% Plot
FS=16;
figure(1)
set(gcf,'Units','normalized','OuterPosition',[0 0 .25 .5])
semilogy(t,N.e,'k-',t,N.O2p,'k--',t,N.CO2p,'k-.',t,N.Op,'k:')
xlabel('t (s)','FontSize',FS);
ylabel('N_\alpha (#)','FontSize',FS);
legend('e','O_2^+','CO_2^+','O^+','Location','best')
legend('Boxoff')
set(gca,'TickDir','Out','FontSize',FS)

figure(2)
set(gcf,'Units','normalized','OuterPosition',[.25 0 .25 .5])
plot(t,(N.O2p+N.CO2p+N.Op+N.e)./(N.O2p(1)+N.CO2p(1)+N.Op(1)+N.e(1)),'k-')
ylim([.95 1])
xlim([0 375])
xlabel('t (s)','FontSize',FS);
ylabel('\Sigma N_\alpha(t) / \Sigma N_\alpha(0)','FontSize',FS);
set(gca,'TickDir','Out','FontSize',FS)
print(2,'-dps','20nTNeu.ps')