R_mer = 57910000000;
M_mer = 3.285*(10^23);
V_mer = 47870;

R_ven = 108200000000;
M_ven = 4.867*(10^24);
V_ven = 35020;

R_tel = 149600000000;
M_tel = 5.972*(10^24);
V_tel = 29783;

R_mar = 227900000000;
M_mar = 6.39*(10^23);
V_mar = 24077;

M_hel = 1.989*(10^30);
V_hel = -1*(V_mer*M_mer + V_ven*M_ven + V_tel*M_tel + V_mar*M_mar)/M_hel;

r_hel = (R_mer*M_mer + R_ven*M_ven + R_tel*M_tel + R_mar*M_mar)...
	/(M_hel - M_mer - M_ven - M_tel - M_mar);



G = 6.67*(10^(-11));
m = [M_hel, M_mer, M_ven, M_tel, M_mar];
x0 = [0, R_mer, R_ven, R_tel, R_mar];
y0 = [0, 0, 0, 0, 0];
vx0 = [0, 0, 0, 0, 0];
vy0 = [V_hel, V_mer, V_ven, V_tel, V_mar];
tmax = 12000000000;
dt = 10000;

[x,y,vx,vy,ax,ay,t] = orbit_Nbody(G,m,x0,y0,vx0,vy0,dt,tmax);

figure(1);
for i = 1:length(m)
   plot(x(:,i),y(:,i),'-o', 'MarkerSize', 1.2);
   hold on;
end
plot(0,0,'-o', 'MarkerSize', 2.2);
axis equal;

[periodtime, periods] = period(x, t)

figure(2);
for i = 1:length(m)
    plot(i, periodtime(i), '-o', 'MarkerSize', 4);
    hold on;
end

[ekt, ept, ett] = energy(G,m,x,y,vx,vy);

figure(3);
plot(t,ett, 'r');
hold on;
plot(t+0.1,ept, 'b');
hold on;
plot(t+0.2,ekt + 5*(10^31), 'g');

[p, ptot] = momentum(m, vx, vy);
figure(4);
plot(t, ptot);

