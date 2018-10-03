clc
close all

listing = dir("/home/alex/Dropbox/CC298/proj01/data/");

% x = load("xgrid_ni21nj12.dat");
% y = load("ygrid_ni21nj12.dat");

x = load("xgrid_ni105nj12.dat");
y = load("ygrid_ni105nj12.dat");

res = load("residual.dat");

p = load("p.dat");
pref1 = csvread("pressure_martins_azevedo.csv");
pref1(:,1) = pref1(:,1) + abs(pref1(1,1));

pref2 = csvread("pressure_maciel.csv");

ft_to_m = 0.3048;
pref2(:,1) = pref2(:,1).*ft_to_m;

nmax = 60;
xref = linspace(0,max(x(:,1)),nmax);

figure
title("Residual")
semilogy(res(:,1),res(:,2))
grid on

figure
title("Wall Pressure")
plot(x(:,1),p(:,1),pref2(:,1),pref2(:,2),"dr")
xlabel("x/l")
ylabel("p/pt")

