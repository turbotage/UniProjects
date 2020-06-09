
%inte länge i omloppsbana någonstans mellan 1.4 och 1.45 ... sqrt(2)?

G = 1;
M = 10;
m = 0.01;
x0 = 10;
y0 = 0;
vx0 = 0;
vy0 = 0.75;
tmax = 300;
dt = 0.001;

%simulate
[x,y,vx,vy,ax,ay,t] = orbit_1body(G,M,x0,y0,vx0,vy0,dt,tmax);

%plot orbits
figure(1);
for i = 1:length(m)
   plot(x(:,i),y(:,i),'-o', 'MarkerSize', 1.2);
   hold on;
   quiver(x(1:50:end,i),y(1:50:end,i),vx(1:50:end,i),vy(1:50:end,i), 0.5);
   hold on;
   quiver(x(1:50:end,i),y(1:50:end,i),ax(1:50:end,i),ay(1:50:end,i), 0.5);
   hold on;
end
plot(0,0,'-o', 'MarkerSize', 10);

mx = (x*m)./(M+m);
my = (y*m)./(M+m);
hold on;
plot(mx,my,'-*');
axis equal;
title('Orbits');
legend('satellite pos','satellite vel','satellite acc','planet','centre of mass');
xlabel('x [m]');
ylabel('y [m]');


%calc periods
[periodtime, periods, max_periods] = period(x,t);
%plot periods
figure(2);
plot([1:max_periods], periodtime,'o', 'MarkerSize', 4); 
xlabel('period');
ylabel('T [s]');
title('Periodtime');

%calc energies
ekt = m*(vx.*vx + vy.*vy)*0.5;
ept = -G*M*m./sqrt(x.*x + y.*y);
ett = ekt + ept;
%plot energies
figure(3);
plot(t,ett, 'r');
hold on;
plot(t,ept, 'b');
hold on;
plot(t,ekt, 'g');
legend('Total Energy','Potential Energy','Kinetic Energy');
xlabel('t [s]');
ylabel('E [J]');
title('Energy');

%plot momentum
[p, ptot] = momentum(m, vx, vy);
figure(4);
plot(t, ptot);
legend('Momentum');
xlabel('t [s]');
ylabel('p [kgm/s]');
title('Momentum');

