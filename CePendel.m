A = importdata("CePendel2.csv");

lengths = A(:,1);
times = A(:,2);

x = lengths;
y = times;

plot(x, y);
print -deps epsFig2

for i = 1:length(lengths)
    y(i) = lengths(i)*times(i)*times(i)/(4*pi*pi); %the y in y = kx + m
    x(i) = lengths(i)*lengths(i); %the x in y = kx + m
end

x = x.';
y = y.';

p = polyfit(x,y,1)
figure(1);
plot(x,y, '-o')
print -deps epsFig

%column 1 is pendulum lengths, column 2 is the angle and column 3 is 
B = importdata("CePendel1Angle.csv");
lengths = B(:,1);
angles = B(:,2);
times = B(:,3);

%times for small angles
tForSmallAngles = 2*pi*sqrt((p(2)/lengths(1))+(p(1)*lengths(1)));

for i = 1:length(times)
   times(i) = times(i)/tForSmallAngles; 
end

p = polyfit(angles, times, 3)
figure(2);
plot(angles, times, '-o');
print -deps epsFig2 
