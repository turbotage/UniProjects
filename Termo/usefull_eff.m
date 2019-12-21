

Q_h = -2.3221e+04*10^(-4);
Q_c = 1.9582e+04*10^(-4);

m = [0.1,0.2,0.3,0.4,0.5];
F = [2.0, 4.0, 5.7, 7.2,9.2];
strob_f = [16.8, 12, 15, 14, 11];
dots = [4,3,4,4,4];
f = strob_f ./ dots
w = 2*pi*f;

r = 0.05;

P_out = r.*w.*(F);
P_in = 11*12.25;



index = find(P_out/P_in == max(P_out/P_in));
W_e = m(index)*9.82;

size1 = (P_in - Q_h*f(index))/P_in
size2 = abs((Q_c*f(index))/P_in)
size3 = abs((W_e*f(index) - P_out(index))/P_in)
size4 = abs(P_out(index)/P_in)

sizeTot = size1 + size2 + size3 + size4
