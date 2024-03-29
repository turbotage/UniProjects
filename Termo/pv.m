p = [0, 10,20, 35, 50, 60, 75, 85, 95, 100.5,115, ...
    100.5 ,95, 85, 75, 60, 50, 35, 20, 4.667, 0];

p = p + 101.325;

V = [270, 230.377, 214.422, 194.453, 177.547, 165.396,...
    156.415, 149.547, 140.566, 132.626, 130, 161.699,166.45, 179.13, 192.34,...
    203.96, 218.23, 234.60, 251.5,270,270];

k = boundary(p.',V.');

plot(V,p,'-b');
title('pV-diagram for hot air engine');
xlabel('V [cm^3]');
ylabel('p [kPa]');
legend('pV-diagram');

T_c = 273.15 + 20;
%hold on;
%plot(V(k),p(k),'-r');

W_all = polyarea(V,p);

W_34 = trapz(V(11:20), p(11:20))
W_12 = trapz(V(1:11),p(1:11))

T_h = -W_34*T_c/W_12 - 273.15

eta = (W_34 + W_12)/W_34;




