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
[topperPath, Folder] = fileparts(pwd);
cd viz_fun

%% Load Cusp data
raw = importdata('../../NewCusp/viz_dir/diagnostics.dat');

%% Resize data
% S.N = size(N.data);
Cusp.t         = raw.data(:,1);
N              = length(Cusp.t);
Cusp.Fi.s(3,N) = 0.0;
Cusp.Fi.n(3,N) = 0.0;
Cusp.Fi.w(3,N) = 0.0;
Cusp.Fi.e(3,N) = 0.0;
Cusp.Fi.b(3,N) = 0.0;
Cusp.Fi.t(3,N) = 0.0;
Cusp.Fe.w(1,N) = 0.0;
Cusp.Fe.e(1,N) = 0.0;
Cusp.Fe.s(1,N) = 0.0;
Cusp.Fe.n(1,N) = 0.0;
Cusp.Fe.b(1,N) = 0.0;
Cusp.Fe.t(1,N) = 0.0;

for l=1:3
    Cusp.Fi.w(l,:) = raw.data(:,3+(l-1)*7);
    Cusp.Fi.e(l,:) = raw.data(:,4+(l-1)*7);
    Cusp.Fi.s(l,:) = raw.data(:,5+(l-1)*7);
    Cusp.Fi.n(l,:) = raw.data(:,6+(l-1)*7);
    Cusp.Fi.b(l,:) = raw.data(:,7+(l-1)*7);
    Cusp.Fi.t(l,:) = raw.data(:,8+(l-1)*7);
end
l=4;
Cusp.Fe.w(1,:) = raw.data(:,3+(l-1)*7);
Cusp.Fe.e(1,:) = raw.data(:,4+(l-1)*7);
Cusp.Fe.s(1,:) = raw.data(:,5+(l-1)*7);
Cusp.Fe.n(1,:) = raw.data(:,6+(l-1)*7);
Cusp.Fe.b(1,:) = raw.data(:,7+(l-1)*7);
Cusp.Fe.t(1,:) = raw.data(:,8+(l-1)*7);

%% Load Arcades data
raw = importdata('../../NewArcades/viz_dir/diagnostics.dat');

%% Resize data
% S.N = size(N.data);
Arca.t         = raw.data(:,1);
N              = length(Arca.t);
Arca.Fi.s(3,N) = 0.0;
Arca.Fi.n(3,N) = 0.0;
Arca.Fi.w(3,N) = 0.0;
Arca.Fi.e(3,N) = 0.0;
Arca.Fi.b(3,N) = 0.0;
Arca.Fi.t(3,N) = 0.0;
Arca.Fe.w(1,N) = 0.0;
Arca.Fe.e(1,N) = 0.0;
Arca.Fe.s(1,N) = 0.0;
Arca.Fe.n(1,N) = 0.0;
Arca.Fe.b(1,N) = 0.0;
Arca.Fe.t(1,N) = 0.0;

for l=1:3
    Arca.Fi.w(l,:) = raw.data(:,3+(l-1)*7);
    Arca.Fi.e(l,:) = raw.data(:,4+(l-1)*7);
    Arca.Fi.s(l,:) = raw.data(:,5+(l-1)*7);
    Arca.Fi.n(l,:) = raw.data(:,6+(l-1)*7);
    Arca.Fi.b(l,:) = raw.data(:,7+(l-1)*7);
    Arca.Fi.t(l,:) = raw.data(:,8+(l-1)*7);
end
l=4;
Arca.Fe.w(1,:) = raw.data(:,3+(l-1)*7);
Arca.Fe.e(1,:) = raw.data(:,4+(l-1)*7);
Arca.Fe.s(1,:) = raw.data(:,5+(l-1)*7);
Arca.Fe.n(1,:) = raw.data(:,6+(l-1)*7);
Arca.Fe.b(1,:) = raw.data(:,7+(l-1)*7);
Arca.Fe.t(1,:) = raw.data(:,8+(l-1)*7);

%% Plot
FS=10;
figure(1)
set(gcf,'Units','normalized','OuterPosition',[.25 0 .25 1],'Color',[1 1 1])

minF = -1e23;
maxF = 5e23;

subplot(2,1,1)
plot(...
    Arca.t,Arca.Fi.t(1,:),'k+-' ,Arca.t,Arca.Fi.t(2,:),'kx-' ,Arca.t,Arca.Fi.t(3,:),'k*-'  ...
    ,...
    Cusp.t,Cusp.Fi.t(1,:),'b+--',Cusp.t,Cusp.Fi.t(2,:),'bx--',Cusp.t,Cusp.Fi.t(3,:),'b*--' ...
    );
ylim([minF maxF])
xlabel('t (s)','FontSize',FS);
ylabel('F_\alpha^{top} (#/s)','FontSize',FS);
legend('O_2^+','CO_2^+','O^+','e','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)
axis tight

subplot(2,1,2)
plot(...
    Arca.t,Arca.Fe.t(1,:),'k' ...
    ,...
    Cusp.t,Cusp.Fe.t(1,:),'b' ... 
    );
ylim([minF maxF])
xlabel('t (s)','FontSize',FS);
ylabel('F_e^{top} (#/s)','FontSize',FS);
legend('Arcades','Cusp','location','best')
legend('boxoff')
set(gca,'TickDir','Out','FontSize',FS)
axis tight

print(1,'-depsc','ComparisonTotalFlowThruTop.eps')