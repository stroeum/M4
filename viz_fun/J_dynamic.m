clear all
close all
beep off
clc

%% Model currents

%% Scan files and determine observation times
folder = '../viz_dir/';
J.x.files = dir([folder,'Jx_*.dat']);
J.y.files = dir([folder,'Jy_*.dat']);
J.z.files = dir([folder,'Jz_*.dat']);
N = length(J.x.files);

%% Creates vectors containing times and current density components
t(N) = 0;
n=1;
t(n) = str2num(J.x.files(n).name(4:end-4));
tmp = importdata([folder,J.x.files(n).name]);
J.x.dat = tmp(:,2);
zobs = 175;
z = 100+tmp(:,1);
kobs = find(floor(z)==zobs);
for n=1:N
    t(n) = str2num(J.x.files(n).name(4:end-4));
    tmp = importdata([folder,J.x.files(n).name]);
    J.x.dat = tmp(:,2);
    J.x.max(n) = max(abs(J.x.dat));
    J.x.obs(n) = J.x.dat(kobs);
    tmp = importdata([folder,J.y.files(n).name]);
    J.y.dat = tmp(:,2);
    J.y.max(n) = max(abs(J.y.dat));
    J.y.obs(n) = J.y.dat(kobs);
    tmp = importdata([folder,J.z.files(n).name]);
    J.z.dat = tmp(:,2);
    J.z.max(n) = max(abs(J.z.dat));
    J.z.obs(n) = J.z.dat(kobs);
    
    J.mag.dat = (J.x.dat.^2 + J.y.dat.^2 + J.z.dat.^2).^.5;
    J.mag.max(n) = max(J.mag.dat);
    J.mag.obs(n) = J.mag.dat(kobs);
end

%% Plots
figure(1)
ft = 16;
set(gcf,'Units','normalized','OuterPosition', [0 0 .25 1]);
plot(t,J.x.max,'k--', t,J.y.max,'k:', t,J.z.max,'k-.', t,J.mag.max,'k-');
xlabel('t (s)','FontSize',ft);
ylabel('Current density (nA/m^2)','FontSize',ft);
set(gca,'FontSize',ft,'TickDir','out','XMinorTick','on','YMinorTick','on');
axis tight
legend('J_x','J_y','J_z','|J|','Location','best')
legend('boxoff')
print(1,'-dps','J_max_bw.ps')

figure(2)
set(gcf,'Units','normalized','OuterPosition', [.25 0 .25 1]);
plot(t,J.x.obs,'k--', t,J.y.obs,'k:', t,J.z.obs,'k-.', t,J.mag.obs,'k-');
xlim([0 125])
xlabel('t (s)','FontSize',ft);
ylabel('Current density (nA/m^2)','FontSize',ft);
set(gca,'FontSize',ft,'TickDir','out','XMinorTick','on','YMinorTick','on');
%axis tight
legend('J_x','J_y','J_z','|J|','Location','best')
legend('boxoff')
print(2,'-dps','J20nT.ps')