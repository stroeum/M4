%% Convert t.dat
close all
clear all
clc

!rm -rfv ../viz_dir/t.dat
s1 = 2.5e3;
s2 = 1e5;

ds = s2/s1;

t_in = load('../output/t.dat');
N = length(t_in);

m=1;
for n=1:N
    if rem(n-1,ds)==0
        t_out(m,1) = t_in(n);
        m = m+1;
    end
end
save('../viz_dir/t.dat','t_out','-ASCII')