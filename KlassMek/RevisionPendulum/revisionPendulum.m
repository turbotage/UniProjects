format long;

M = load('iter3.csv'); % Meassurements

L = M(:,1) .* 10^(-3);
T_f = M(:,2) .* 10^(-6);
T_r = M(:,3) .* 10^(-6);

%approx 1 275 mm

%0.9535

[P_f, S_f] = polyfit(L,T_f,1);
[P_r, S_r] = polyfit(L,T_r,1);

[T_fit_f, delta_f] = polyval(P_f, L, S_f);
[T_fit_r, delta_r] = polyval(P_r, L, S_r);

%l1 = (P_f(2) - 2*delta_f) - (P_r(2) + 2*delta_r) / (P_r(1) - P_f(1));
%l2 = (P_f(2) + 2*delta_f) - (P_r(2) - 2*delta_r) / (P_r(1) - P_f(1));
%disp(l1)
%disp(l2)

figure(1);
%T_f
plot(L,T_fit_f + 2*delta_f,'b', L, T_fit_f - 2*delta_f, 'b');
hold on;
plot(L, T_fit_f);
hold on;
plot(L, T_f, '.');
hold on;
%T_r
plot(L,T_fit_r + 2*delta_r,'r',L, T_fit_r - 2*delta_r, 'r');
hold on;
plot(L, T_fit_r);
hold on;
plot(L, T_r, '.');


%plot(L,T_f); % T_f
%hold on;
%plot(L,T_r); % T_r

l_min = 0.95121436;
l_max = 0.9546783;
l_mid = 0.95296761; 

t_min = 1.956236606;
t_max = 1.9580200287;
t_mid = 1.95713928;

l_matrix = [l_min, l_max, l_mid];
t_matrix = [t_min, t_max, t_mid];

g = l_matrix .* (2*pi./t_matrix).^2

g_mean = (g(1)+g(2))/2
g_uncertainty = (g(2)-g(1))/2