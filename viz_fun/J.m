clear all
close all
beep off
clc

%% Model currents
tmp = load('Jx.dat');
z = tmp(:,1)+100;
tmp = load('Jx.dat');
J.x = tmp(:,2);
tmp = load('Jy.dat');
J.y = tmp(:,2);
tmp = load('Jz.dat');
J.z = tmp(:,2);
J.o = sqrt(J.x.^2+J.y.^2+J.z.^2);

set(gcf,'Units','normalize','Color','white','OuterPosition',[.25 .25 .25 .75])
plot(J.x,z,'r',J.y,z,'b',J.z,z,'g',J.o,z,'k')
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',18);
xlabel('Current density (nA/m^2)')
ylabel('z (km)')
legend('J_x','J_y','J_z','|J|','Location','best')
legend('boxoff')

%% Theoretical currents
%%
b.u    = 20e-9; % T
%b.u    = 2000e-9; % T

set(gcf,'Units','normalize','Color','white')
q     = 1.60217646e-19; % C
m.O2p = 5.3135e-26; % kg
m.CO2 = 7.3079e-26; % kg
m.e   = 9.3109e-31; % kg

tmp = load('ne.dat');
z   = tmp(:,1)+100; %_km
nz  = length(z);
tmp = load('ne.dat');
n.e = tmp(:,2); %_m^-3

tmp = load('Te.dat');
T.e = tmp(:,2); %_K
tmp = load('TO2p.dat');
T.i = tmp(:,2); %_K

v     = importdata('../data/MGS/Ls180_LT14_MY24_solarmod.dat');
v.data(any(isnan(v.data),2),:) =[];
s.z   = v.data(:,1);     %_km
s.n.g = v.data(:,2)*1e6; %_m-3
s.T.g = v.data(:,9);     %_K
n.g   = interp1(s.z,s.n.g,z); %_m-3
T.g   = interp1(s.z,s.T.g  ,z); %_K

%% Cyclotron frequencies
Omega.i = q*b.u/(2*pi*m.O2p);
Omega.e = q*b.u/(2*pi*m.e);

%% Collision frequencies
for k=1:nz
    v.in(k)  = Collision_Frequency('O2+','CO2' ,n.e(k), n.g(k), T.i(k), T.g(k), z(k));
    v.en(k)  = Collision_Frequency('e'  ,'CO2' ,n.e(k), n.g(k), T.e(k), T.g(k), z(k));
end

k1=min(find(v.en<=Omega.e));
k2=max(find(v.in>=Omega.i));

for k=1:nz
    if(z(k)>z(k2))
        V.e.x(k) = 0;
        V.e.y(k) = 100;
        V.i.x(k) = 0;
        V.i.y(k) = 100;
        V.n(k) = 100;
    elseif(z(k)>=z(k1) && z(k)<=z(k2))
        V.e.x(k) = 0;
        V.e.y(k) = 100;
        V.i.x(k) = 100;
        V.i.y(k) = 0;
        V.n(k) = 100;
    else
        V.e.x(k) = 100;
        V.e.y(k) = 0;
        V.i.x(k) = 100;
        V.i.y(k) = 0;
        V.n(k) = 100;
    end
end

figure(2)
subplot(121)
plot(V.e.x,z,'r',V.e.y,z,'b')
subplot(122)
plot(V.i.x,z,'r',V.i.y,z,'b')

figure(3)
th.J.x = -q*n.e'.*(v.en./Omega.e.*V.e.y       + v.in./Omega.i.*V.i.y); 
th.J.y = -q*n.e'.*(v.en./Omega.e.*(V.n-V.e.x) + v.in./Omega.i.*(V.n-V.i.x));
% plot(th.J.x*1e9,z,'r',th.J.y*1e9,z,'b+')

set(gcf,'Units','normalize','Color','white','OuterPosition',[.25 .25 .25 .75])
plot(J.x,z,'r',J.y,z,'b',J.z,z,'g',J.o,z,'k',th.J.x*1e9,z,'r+',th.J.y*1e9,z,'b+')
set(gca,'XminorTick','on','YMinorTick','on','Tickdir','out','FontSize',18);
xlabel('Current density (nA/m^2)')
ylabel('z (km)')
legend('J_x','J_y','J_z','|J|','Location','best')
legend('boxoff')