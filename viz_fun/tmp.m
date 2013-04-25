clear all
close all
clc
%% Load data
DiR.N.data = load('../../Ext20nTDiRichlet/viz_dir/CarriersCounts.dat');
DiR.t = load('../../Ext20nTDiRichlet/output/t.dat');
Neu.N.data = load('../../20nT/viz_dir/CarriersCounts.dat');
Neu.t = load('../../20nT/output/t.dat');

%% Resize data
S.N = size(DiR.N.data);
DiR.N.e    = DiR.N.data(:,1);
DiR.N.O2p  = DiR.N.data(:,2);
DiR.N.CO2p = DiR.N.data(:,3);
DiR.N.Op   = DiR.N.data(:,4);
S.t = length(DiR.t);
if (S.N(1)<S.t)
    DiR.t = DiR.t(1:S.N(1));
else
    DiR.N.e    = DiR.N.e(1:S.t);
    DiR.N.O2p  = DiR.N.O2p(1:S.t);
    DiR.N.CO2p = DiR.N.CO2p(1:S.t);
    DiR.N.Op   = DiR.N.Op(1:S.t);
end

S.N = size(Neu.N.data);
Neu.N.e    = Neu.N.data(:,1);
Neu.N.O2p  = Neu.N.data(:,2);
Neu.N.CO2p = Neu.N.data(:,3);
Neu.N.Op   = Neu.N.data(:,4);
S.t = length(Neu.t);
if (S.N(1)<S.t)
    Neu.t = Neu.t(1:S.N(1));
else
    Neu.N.e    = Neu.N.e(1:S.t);
    Neu.N.O2p  = Neu.N.O2p(1:S.t);
    Neu.N.CO2p = Neu.N.CO2p(1:S.t);
    Neu.N.Op   = Neu.N.Op(1:S.t);
end


%% Plot
FS=16;
figure(1)
set(gcf,'Units','normalized','OuterPosition',[0 0 .25 .5])
semilogy(DiR.t,DiR.N.e,'k-',DiR.t,DiR.N.O2p,'k--',DiR.t,DiR.N.CO2p,'k-.',DiR.t,DiR.N.Op,'k:')
xlabel('t (s)','FontSize',FS);
ylabel('N_\alpha (#)','FontSize',FS);
legend('e','O_2^+','CO_2^+','O^+','Location','best')
legend('Boxoff')
set(gca,'TickDir','Out','FontSize',FS)

figure(2)
set(gcf,'Units','normalized','OuterPosition',[.25 0 .25 .5])
plot(DiR.t,(DiR.N.O2p+DiR.N.CO2p+DiR.N.Op+DiR.N.e)./(DiR.N.O2p(1)+DiR.N.CO2p(1)+DiR.N.Op(1)+DiR.N.e(1)),'k-',...
     Neu.t,(Neu.N.O2p+Neu.N.CO2p+Neu.N.Op+Neu.N.e)./(Neu.N.O2p(1)+Neu.N.CO2p(1)+Neu.N.Op(1)+Neu.N.e(1)),'k--')
ylim([.95 1])
xlim([0 375])
xlabel('t (s)','FontSize',FS);
ylabel('\Sigma N_\alpha(t) / \Sigma N_\alpha(0)','FontSize',FS);
set(gca,'TickDir','Out','FontSize',FS)
legend('DiRichlet','Neumann');
legend('boxoff')
print(2,'-dps','20nT.ps')