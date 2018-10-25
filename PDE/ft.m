
%% assignment 1
y = [8,12,-4,4];
z = mydft(y)

%% assignment 2 (b)

f = @(x)3 + 2.*cos(45.*x) + 4.*sin(120.*x);
N = 2^8;
x = zeros(1,N);
% (-pi,pi)
for i=0:(N-1)
    x(i+1) = (2*pi*i/N);
end

y = f(x)
z = mydft(y)


figure(1);
subplot(2,2,1)
stem([0:N-1],real(z));
title('Subplot 1: real(z)')

subplot(2,2,2)
stem([0:N-1],imag(z));
title('Subplot 2: imag(z)')

subplot(2,2,3)
stem([0:N-1],abs(z));
title('Subplot 3: abs(z)')

subplot(2,2,4)
plot(x,y);
title('Subplot 4: f(x)')

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
N = 2^6;
x = zeros(1,N);
% (0,2pi)
for i=0:(N-1)
    x(i+1) = (2*pi*i/N);
end

y = g(x)
z = mydft(y)
y_t = myifft(z,length(z))

figure(1);
subplot(2,2,1)
stem([0:N-1],real(z));
title('Subplot 1: real(z)')


subplot(2,2,2)
stem([0:N-1],imag(z));
title('Subplot 2: imag(z)')

subplot(2,2,3)
stem([0:N-1],abs(z));
title('Subplot 3: abs(z)')

subplot(2,2,4)
plot(x,y);
title('Subplot 4: g(x)')

%% assignment 3
y = [8,12,-4,4];
z1 = myfft(y,4)
z2 = mydft(y);
z3 = myifft(z1,4)

