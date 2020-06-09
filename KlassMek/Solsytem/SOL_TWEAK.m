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

R_lun = R_tel - 384400000;
M_lun = M_tel*0.0123;
V_lun = V_tel - 1022;


R_hel = 0;
M_hel = 1.989*(10^30);
V_hel = -1*(V_mer*M_mer + V_ven*M_ven + V_tel*M_tel + V_mar*M_mar + ...
    V_lun*M_lun)/M_hel;

r_offset = (M_mer*R_mer + M_ven*R_ven + M_tel*R_tel + M_mar*R_mar + ...
    M_lun*R_lun)/(M_hel + M_mar + M_ven + M_tel + M_mar + M_lun);




G = 6.67*(10^(-11));
m = [M_hel, M_mer, M_ven, M_tel, M_mar, M_lun];
x0 = [R_hel, R_mer, R_ven, R_tel, R_mar, R_lun] - r_offset;
y0 = [0, 0, 0, 0, 0, 0];
vx0 = [0, 0, 0, 0, 0, 0];
vy0 = [V_hel, V_mer, V_ven, V_tel, V_mar, V_lun];
tmax = 100000000;
dt = 300;

[x,y,vx,vy,ax,ay,t] = orbit_Nbody(G,m,x0,y0,vx0,vy0,dt,tmax);

%%
figure(1);
for i = 1:length(m)
   plot(x(:,i),y(:,i),'-o', 'MarkerSize', 1.2);
   hold on;
end
plot(0,0,'-o', 'MarkerSize', 2.2);
axis equal;

%%
%calculate periodtime and number of periods
[periodtime, periods, max_periods] = period(x, t);
[periodtime_temp, periods_temp, max_periods_temp] = period(x(:,end)-x(:,end-2),t);
periodtime(:,end) = periodtime_temp(1:max_periods);
periods(:,end) = periods_temp;
%plot periodtimes
figure(2);
for i = 1:length(m)
    plot([1:max_periods], periodtime(:,i), 'o', 'MarkerSize', 4);
    hold on;
end
title('Periodtime');

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

disp(periodtime(1,:)./(3600*24));

% 159.6354   88.6424  225.4028  365.5938  683.2535   27.3056
%            88       225       365       687        27

