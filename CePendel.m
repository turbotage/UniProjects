A = importdata("CePendel1Angle.txt");
B = A(:,1);
C = A(:,2);

B = B.';
C = C.';

p = polyfit(B,C,2)

plot(B,C)