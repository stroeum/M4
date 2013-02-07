clear all
close all
clc

%% V1 measurements
dir_path='../data/';

ft = fittype( 'a*exp(1-(x-b)/c-exp(-(x-b)/c))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display    = 'Off';

z = 100:1:400;

%% O2+
V1.O2p    = load([dir_path,'Viking','/O2p' ,'.dat']);
V1.n.O2p  = V1.O2p(:,1) *1e6;
V1.z.O2p  = V1.O2p(:,2) ;

zmin = 120;
zmax = 390;
xdata = V1.z.O2p(zmin<V1.z.O2p<zmax);
ydata = V1.n.O2p(zmin<V1.z.O2p<zmax);

opts.Lower      = [1e9 100 1];
opts.StartPoint = [6e10 120 10];
opts.Upper      = [1e12 300 150];
[x, gof]        = fit( xdata, ydata, ft, opts );

disp(['x.a.O2p  = ',num2str(x.a,'%8.3e'),';'])
disp(['x.b.O2p  = ',num2str(x.b,'%8.3e'),';'])
disp(['x.c.O2p  = ',num2str(x.c,'%8.3e'),';'])

A.O2p = x.a;
Z.O2p = (z-x.b)/x.c;
O2p = A.O2p*exp(1-Z.O2p-exp(-Z.O2p));

%% CO2+
V1.CO2p   = load([dir_path,'Viking','/CO2p','.dat']);
V1.n.CO2p = V1.CO2p(:,1)*1e6;
V1.z.CO2p = V1.CO2p(:,2);

zmin = 100;
zmax = 250;
xdata = V1.z.CO2p(zmin<V1.z.CO2p<zmax);
ydata = V1.n.CO2p(zmin<V1.z.CO2p<zmax);

opts.Lower      = [1e9 100 1];
opts.StartPoint = [.75e10 150 10];
opts.Upper      = [1e11 300 150];
[x, gof]        = fit( xdata, ydata, ft, opts );
c(1) = x.a;
c(2) = x.b;
c(3) = x.c;
disp(['x.a.CO2p = ',num2str(x.a,'%8.3e'),';'])
disp(['x.b.CO2p = ',num2str(x.b,'%8.3e'),';'])
disp(['x.c.CO2p = ',num2str(x.c,'%8.3e'),';'])

A.CO2p = x.a;
Z.CO2p = (z-x.b)/x.c;
CO2p = A.CO2p*exp(1-Z.CO2p-exp(-Z.CO2p));

%% O+
V1.Op     = load([dir_path,'Viking','/Op'  ,'.dat']);
V1.n.Op   = V1.Op(:,1)  *1e6;
V1.z.Op   = V1.Op(:,2)  ;

zmin = 180;
zmax = 280;
xdata = V1.z.Op(zmin<V1.z.Op<zmax);
ydata = V1.n.Op(zmin<V1.z.Op<zmax);

opts.Lower      = [1e8 200 1];
opts.StartPoint = [5e8 220 10];
opts.Upper      = [1e9 300 25];
[x, gof]        = fit( xdata, ydata, ft, opts );

disp(['x.a.Op   = ',num2str(x.a,'%8.3e'),';'])
disp(['x.b.Op   = ',num2str(x.b,'%8.3e'),';'])
disp(['x.c.Op   = ',num2str(x.c,'%8.3e'),';'])

A.Op = x.a;
Z.Op = (z-x.b)/x.c;
Op   = A.Op*exp(1-Z.Op-exp(-Z.Op));

set(gcf,'Units','Normalized','OuterPosition',[0 0 .25 .5])
semilogx(O2p,z,'b-', V1.n.O2p, V1.z.O2p, 'bo', CO2p,z,'r-', V1.n.CO2p, V1.z.CO2p, 'r^', Op, z, 'g-', V1.n.Op, V1.z.Op, 'gs')
xlabel('n_\alpha (cm^{-3}')
ylabel('z (km)')
xlim([1e5 1e11])
ylim([100 400])

