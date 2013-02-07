clear all
close all
clc

eps = 1e-10;

%%
% b.u    = 20e-9; % T
b.u    = 2000e-9; % T

set(gcf,'Units','normalize','Color','white')
gamma  = 5/3;
c      = 299792458;  % m/s
kB     = 1.3806488e-23; % J/K
mu.o   = 4*pi*1e-7; % V.s/(A.m)
q      = 1.60217646e-19; % C
m.O2p  = 5.3135e-26; % kg
m.CO2p = 7.3079e-26; % kg
m.Op   = 2.6567e-26; % kg
m.CO2  = 7.3079e-26; % kg
m.O    = 2.6567e-26; % kg
m.N2   = 4.6833e-26; % kg
m.CO   = 4.6833e-26; % kg
m.e    = 9.3109e-31; % kg

folder = '../viz_dir/';
tmp    = load([folder,'ne.dat']);
z      = tmp(:,1)+100; %_km
nz     = length(z);
tmp    = load([folder,'ne.dat']);
n.e    = tmp(:,2); %_m^-3
tmp    = load([folder,'nO2p.dat']);
n.O2p  = tmp(:,2); %_m^-3
tmp    = load([folder,'nCO2p.dat']);
n.CO2p = tmp(:,2); %_m^-3
tmp    = load([folder,'nOp.dat']);
n.Op   = tmp(:,2); %_m^-3

tmp    = load([folder,'pe.dat']);
p.e    = tmp(:,2); %_Pa
tmp    = load([folder,'pO2p.dat']);
p.O2p  = tmp(:,2); %_Pa
tmp    = load([folder,'pCO2p.dat']);
p.CO2p = tmp(:,2); %_Pa
tmp    = load([folder,'pOp.dat']);
p.Op   = tmp(:,2); %_Pa

tmp    = load([folder,'Te.dat']);
T.e    = tmp(:,2); %_K
tmp    = load([folder,'TO2p.dat']);
T.O2p  = tmp(:,2); %_K
tmp    = load([folder,'TCO2p.dat']);
T.CO2p = tmp(:,2); %_K
tmp    = load([folder,'TOp.dat']);
T.Op   = tmp(:,2); %_K

v       = importdata('../../../Current Work/Profiles/data/MGS/Ls180_LT14_MY24_solarmod.dat');
v.data(any(isnan(v.data),2),:) =[];
s.z     = v.data(:,1);     %_km
s.n.CO2 = v.data(:,2)*1e6; %_m-3
s.n.O   = v.data(:,3)*1e6; %_m-3
s.n.N2  = v.data(:,4)*1e6; %_m-3
s.n.CO  = v.data(:,5)*1e6; %_m-3
s.T.g   = v.data(:,9);     %_K

n.CO2   = interp1(s.z,s.n.CO2,z); %_m-3
n.O     = interp1(s.z,s.n.O  ,z); %_m-3
n.N2    = interp1(s.z,s.n.N2 ,z); %_m-3
n.CO    = interp1(s.z,s.n.CO ,z); %_m-3
T.g     = interp1(s.z,s.T.g  ,z); %_K

beta    = beta_eff(T.e);   %_cm6/s

n.g     = n.CO2 + n.O + n.N2 + n.CO; %_m^-3
rho.g   = m.CO2*n.CO2 + m.O*n.O + m.N2*n.N2 + m.CO*n.CO;
rho.c   = m.e*n.e + m.O2p*n.O2p + m.CO2p*n.CO2p + m.Op*n.Op;
p.g     = n.g*kB.*T.g;

T.r     = (T.g + T.O2p)/2;

%% Cyclotron frequencies
Omega(1) = q*b.u/(2*pi*m.O2p);
Omega(2) = q*b.u/(2*pi*m.CO2p);
Omega(3) = q*b.u/(2*pi*m.Op);
Omega(4) = q*b.u/(2*pi*m.e);


%% Collision frequencies
vcn  = zeros(10,nz);
for k=1:nz
    vcn(1,k)  = Collision_Frequency('O2+' ,'CO2' ,n.O2p(k) ,n.CO2(k) ,T.O2p(k) ,T.g(k)   ,z(k));
    vcn(2,k)  = Collision_Frequency('CO2+','CO2' ,n.CO2p(k),n.CO2(k) ,T.CO2p(k),T.g(k)   ,z(k));
    vcn(3,k)  = Collision_Frequency('e'   ,'CO2' ,n.e(k)   ,n.CO2(k) ,T.e(k)   ,T.g(k)   ,z(k));
    vcn(4,k)  = Collision_Frequency('e'   ,'e'   ,n.e(k)   ,n.e(k)   ,T.e(k)   ,T.e(k)   ,z(k));
    vcn(5,k)  = Collision_Frequency('e'   ,'O2+' ,n.e(k)   ,n.O2p(k) ,T.e(k)   ,T.O2p(k) ,z(k));
    vcn(6,k)  = Collision_Frequency('e'   ,'CO2+',n.e(k)   ,n.CO2p(k),T.e(k)   ,T.CO2p(k),z(k));
    vcn(7,k)  = Collision_Frequency('O2+' ,'O2+' ,n.O2p(k) ,n.O2p(k) ,T.O2p(k) ,T.Op(k)  ,z(k));
    vcn(8,k)  = Collision_Frequency('CO2+','O2+' ,n.CO2p(k),n.O2p(k) ,T.CO2p(k),T.O2p(k) ,z(k));
    vcn(9,k)  = Collision_Frequency('O2+' ,'CO2+',n.O2p(k) ,n.CO2p(k),T.Op(k)  ,T.CO2p(k),z(k));
    vcn(10,k) = Collision_Frequency('CO2+','CO2+',n.CO2p(k),n.CO2p(k),T.CO2p(k),T.CO2p(k),z(k));
end

%% Sound and Alfven velocities
v.s = sqrt(gamma*p.g./rho.g);
v.A = b.u./sqrt(mu.o*rho.c);

%% Plots
figure(1)
set(gcf,'Units','normalized','Color',[1 1 1],'OuterPosition',[.5,0,.25,1])
semilogx(...
    vcn(1,:)    ,z,'r+-' ,vcn(2,:)   ,z,'gp-',vcn(3,:)    ,z,'b^-' ,...
    vcn(4,:)    ,z,'rh-' ,vcn(5,:)   ,z,'g-' ,vcn(6,:)    ,z,'bs-',...
    vcn(7,:)    ,z,'r-'  ,vcn(8,:)   ,z,'gv-',vcn(9,:)    ,z,'b-' ,vcn(10,:)   ,z,'k-');
xlim([min(min(vcn)),max(max(vcn))])
ylim([100 max(z)])
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',18);
xlabel('\nu_{\alpha-n} (s^{-1})')
ylabel('z (km)')
xlim([1e-10 1e6])
legend('O_2^+-CO_2','CO_2^+-CO_2','e-CO_2','e-e','e-O_2^+','e-CO_2^+','O_2^+-O_2^+','CO_2^+-O_2^+','O_2^+-CO_2^+','CO_2^+-CO_2^+','Location','best')
legend('boxoff')

figure(2)
set(gcf,'Units','normalized','Color',[1 1 1],'OuterPosition',[.75,0,.25,.75])
subplot(221)
semilogx(n.O2p*1e-6, z, 'b', n.CO2p*1e-6, z, 'r', n.Op*1e-6, z, 'g', n.e*1e-6, z, 'k');
xlim([1e1 2e5])
ylim([100 max(z)])
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',18);
set(gca,'XTick',[1e1 1e2 1e3 1e4 1e5])
xlabel('n_\alpha (cm^{-3})')
ylabel('z (km)')
legend('O_2^+','CO_2^+','O^+','e','Location','best')
legend('boxoff')

subplot(222)
semilogx(n.CO2*1e-6, z, 'r',n.O*1e-6, z, 'b');
xlim([1e4 2e12])
ylim([100 max(z)])
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',18);
set(gca,'XTick',[1e4 1e6 1e8 1e10 1e12])
xlabel('n_\alpha (cm^{-3})')
ylabel('z (km)')
legend('CO_2','O','Location','best')
legend('boxoff')

subplot(223)
semilogx(T.g, z, 'k', T.O2p, z, 'r', T.e, z, 'b', T.r, z, 'k--',[235 235],[z(1) z(end)], 'k:');
%T.O2p*1e-6, z, 'b', T.CO2p*1e-6, z, 'r', T.Op*1e-6, z, 'g', T.e*1e-6, z, 'k');
%xlim([1e1 1e4])
ylim([100 max(z)])
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',18);
xlabel('T_\alpha (K)')
ylabel('z (km)')
legend('neutral','ion','electron','Location','best')
legend('boxoff')

subplot(224)
semilogx(p.e, z, 'k', p.O2p, z, 'b', p.CO2p, z, 'r', p.Op, z, 'g');
xlim([1e-12 1e-9])
ylim([100 max(z)])
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',18);
set(gca,'XTick',[1e-12 1e-11 1e-10 1e-9 1e-8])
xlabel('P_\alpha (Pa)')
ylabel('z (km)')
legend('e','O_2^+','CO_2^+','O^+','Location','best')
legend('boxoff')

%%
figure(3)
k1=min(find(vcn(3,:)<=Omega(4)));
k2=max(find(vcn(1,:)>=Omega(1)));
set(gcf,'Units','normalized','Color',[1 1 1],'OuterPosition',[0,0.1,.33,.7])

subplot('Position',[.1 .15 .5 .7])
xx = [1e-1 1e5];
% xx = [1e-3 1e3];
yy = [100 400];
hold on
rectangle('Position', [xx(1) z(k1) xx(2) z(k2)-z(k1)],'LineStyle','-','EdgeColor','g','FaceColor','g')
semilogx(...
    vcn(1,:)    ,z,'r-' ,vcn(3,:)    ,z,'b-', ...
    Omega(1)*z./z,z,'r-' ,Omega(4)*z./z,z,'b-' );
set(gca,'XScale','log')
hold off
xlim(xx)
ylim(yy)
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',18);
xlabel('\nu_{\alpha-n} (Hz)')
ylabel('z (km)')
legend('O_2^+-CO_2','e-CO_2','\Omega_{O_2^+}','\Omega_{e}','Location','best')
set(gca,'XTick',[1e-1 1e0 1e1 1e2 1e3 1e4 1e5])
legend('boxoff')
box on
title(['B = ',num2str(b.u*1e9),' nT i_z'])

subplot('Position',[.6 .15 .15 .7])
semilogx(n.e*1e-6,z,'k');
xlim([1e3 2e5])
ylim(yy)
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',12,'XAxisLocation','top','FontSize',18);
set(gca,'XTick',[1e3 1e4 1e5])
xlabel('n_e (cm^{-3})')
set(gca,'YTickLabel','');

figure(4)
set(gcf,'Units','normalize','Color','white','OuterPosition',[.5 .25 .25 .75])
semilogx(v.A,z,'k',v.s,z,'k--')
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',18);
ylim([min(z) max(z)]);
xlim([min([min(v.A),min(v.s)]) min(max(max(v.A),max(v.s)),c)]);
xlabel('v_A,c_s (m/s)')
ylabel('z (km)')
legend('v_A','v_s')
legend('boxoff')
title(['B = ',num2str(b.u*1e9),' nT i_z'])

figure(5)
set(gcf,'Units','normalize','Color','white','OuterPosition',[.25 .25 .25 .75])
plot(central_diff(p.O2p,z), z, 'b', central_diff(p.CO2p,z), z, 'r', central_diff(p.Op,z), z, 'g', central_diff(p.e,z), z, 'k');
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',18);
ylim([min(z) max(z)]);
%xlim([min([min(v.A),min(v.s)]) min(max(max(v.A),max(v.s)),c)]);
xlabel('\nabla p_\alpha (Pa/m)')
ylabel('z (km)')
legend('O_2^+','CO_2^+','O^+','e','Location','best')
legend('boxoff')
%title(['B = ',num2str(b.u*1e9),' nT i_z'])

print(2,'-dps','Profiles.ps')
print(3,'-dps','Dynamo.ps')
print(5,'-dps','GradP.ps')
