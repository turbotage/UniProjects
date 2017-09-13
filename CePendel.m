A = importdata("CePendel1Angle.dat");
x = A(:,1);
y = B(:,2);

x = x.';
y = y.';

p = polyfit(x,y,1)
plot(x,y)