close all;
clear all;
clc;

cd ..
[topperPath, Folder] = fileparts(pwd);
cd viz_fun


%% Load data
raw = importdata('../viz_dir/F_CO2p_north_0.dat');
Nz = size(raw,1);
Z = raw(:,1);

F_CO2p_north(1,1:Nz) = raw(:,2);

for n=1:26
    step = (n-1)*10000;
    raw = importdata(['../viz_dir/F_CO2p_north_',num2str(step),'.dat']);
    F_CO2p_north(n+1,1:Nz) = raw(1:Nz,2);
end

FS = 16;
set(gcf,'Units','normalized','OuterPosition',[0 0 .25 1],'Color',[1 1 1])
hold on
for n=1:26
    step = (n-1)*10000; 
    plot(F_CO2p_north(n,:),Z*1e-3,'b-')
    xlim([min(min(F_CO2p_north)) max(max(F_CO2p_north))])
    pause(.1);
end
hold off
xlabel('F_{CO_2^+}^{north}(z) (#)','FontSize',FS);
ylabel('z (km)','FontSize',FS);
legend('Boxoff')
axis tight
set(gca,'TickDir','Out','FontSize',FS)
title([Folder,'@ step ',num2str(step)])
