
%% assignment 1
y_1 = [8,12,-4,4];
z = mydft(y_1);
display(z)
new_y = myidft(z);
display(new_y)

%% assignment 2 (b)

f = @(x)3 + 2.*cos(45.*x) + 4.*sin(120.*x);
N = 2^4;
x = zeros(1,N);
% (0,2pi)
for i=0:(N-1)
    x(i+1) = (2*pi*i/N);
end

y = f(x)
z = mydft(y)


figure(1);
stem([0:N-1],real(z));
title('Subplot 1: real(z)')
xlabel('n')
ylabel('Re z(n)')

figure(2);
stem([0:N-1],imag(z));
title('Subplot 2: imag(z)')
xlabel('n')
ylabel('Im z(n)')

figure(3);
stem([0:N-1],abs(z));
title('Subplot 3: abs(z)')
xlabel('n')
ylabel('Abs z(n)')

figure(4);

plot(x,y);
title('Subplot 4: f(x)')
xlabel('x')
ylabel('f(x)')

a_n = getan(z);
b_n = getbn(z);

%
figure(2);
plot([0:length(a_n)-1],a_n);
title('a_n');

figure(3);
plot([0:((N/2)-2)],b_n);
title('b_n');

%% assignment 2 (c)


g = @(x) abs(cos(x));
N = 2^8;
x = zeros(1,N);

a_exact = zeros(N/2);
a_exact(1) = 4/pi;

for i=1:N/4-1
    a_exact(2*i+1) = -(4/pi*(-1)^i)/(4*i^2-1); 
end

% (0,2pi)
for i=0:(N-1)
    x(i+1) = (2*pi*i/N);
end

y = g(x);
z = mydft(y);
%y_t = myifft(z,length(z))
a_n = getan(z);

figure(1)
stem([0:N-1],real(z))
title('Subplot 1: real(z)')


figure(2)
stem([0:N/2-1],a_n(1:N/2),'*')
hold on
stem([0:N/2-1],a_exact(1:N/2),'or')
hold off
axis([0,120,-0.15,1.3])
xlabel('n')
ylabel('a_n')
legend('Computed a_n','Exact a_n')
title('Plot of computed and exact values of a_n')

figure(3)
stem([0:N-1],imag(z));
title('Subplot 2: imag(z)')

figure(4)
plot(x,y);
title('Subplot 4: g(x)')

%% assignment 2 (d)
h = @(x) x.*((-pi<x) & (x < pi)) + (x - 2*pi).*((pi < x) & (x < 2*pi));

N = 2^4;
x = zeros(1,N);
% (0,2pi)
for i=0:(N-1)
    x(i+1) = (2*pi*i/N);
end
y = h(x);
z=mydft(y);
b_n = getbn(z);
a_n = getan(z);

x2 = 0:0.001:2*pi;
x_i=part_sum(z,x2);


figure(1)
stem([0:N-1],real(z));
hold on;
stem([0:N-1],imag(z));

figure(2)
plot(x2,h(x2), '-b');
hold on;
plot(x2,x_i,'-r');


%% assignment 3
y = [8,12,-4,4];
z1 = myfft(y,4)
z2 = mydft(y);
z3 = myifft(z1,4)

