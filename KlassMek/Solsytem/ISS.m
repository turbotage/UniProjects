%earth
R_tel = 0;
M_tel = 5.972*(10^24);
V_tel = 0;

%satalite
R_sat = 6371000+417500;
M_sat = 419725;
V_sat = 27600/3.6;

%set earth velocity so that V of centre of mass is zero
V_tel = -1*(V_sat*M_sat)/M_tel;

G = 6.67*(10^(-11));
m = [M_tel, M_sat];
x0 = [R_tel, R_sat];
y0 = [0, 0];
vx0 = [0, 0];
vy0 = [V_tel, V_sat];
tmax = 24000;
dt = 0.1;

%simulate
[x,y,vx,vy,ax,ay,t] = orbit_Nbody(G,m,x0,y0,vx0,vy0,dt,tmax);

%plot orbits
figure(1);
for i = 1:length(m)
   plot(x(:,i),y(:,i),'-o', 'MarkerSize', 1.2);
   hold on;
end
plot(0,0,'-o', 'MarkerSize', 10);
axis equal;
title('Orbits');

%calculate periodtime and number of periods
[periodtime, periods, max_periods] = period(x, t);
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
plot(t,ept, 'b');
hold on;
plot(t,ekt, 'g');
title('Energy');

[p, ptot] = momentum(m, vx, vy);
figure(4);
plot(t, ptot);
title('Momentum');

