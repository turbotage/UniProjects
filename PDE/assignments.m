%% assignment 1
clear all;
clc;

y_1 = [8,12,-4,4];
z = mydft(y_1);
display(z)
new_y = myidft(z);
display(new_y)

%% assignment 2 (a)
clear all;
clc;

f = @(x)3 + 2.*cos(15.*x) + 4.*sin(120.*x);
N = 2^8;
x = zeros(1,N);
% (0,2pi)
for i=0:(N-1)
    x(i+1) = (2*pi*i/N);
end

y = f(x);
z = mydft(y);
[a0,a_n,b_n] = myfouriercoeff(z);

figure(1);
plot(x,y);
title('Subplot 4: f(x)')
xlabel('x')
ylabel('f(x)')

figure(2);
plot([0:length(a_n)],[a0,a_n], 'O');
title('a_n');

figure(3);
plot([1:length(b_n)],b_n, 'O');
title('b_n');

% figure(4);
% stem([0:N-1],real(z));
% title('Subplot 1: real(z)')
% xlabel('n')
% ylabel('Re z(n)')
% 
% figure(5);
% stem([0:N-1],imag(z));
% title('Subplot 2: imag(z)')
% xlabel('n')
% ylabel('Im z(n)')
% 
% figure(6);
% stem([0:N-1],abs(z));
% title('Subplot 3: abs(z)')
% xlabel('n')
% ylabel('Abs z(n)')

%% assignment 2 (b)
clear all;
clc;

g = @(x) abs(cos(x));
N = 2^8;
x = zeros(1,N);

a_exact = zeros(N/2);
a_exact(1) = 4/pi;

for i=1:N/4-1
    a_exact(2*i+1) = -(4/pi)*((-1)^i)/(4*i^2-1); 
end

% [0,2pi]
for i=0:(N-1)
    x(i+1) = (2*pi*i/N);
end

y = g(x);
z = mydft(y);
%y_t = myifft(z,length(z))
a_n = getan(z);


figure(1)
stem([0:N/2-1],a_n(1:N/2),'*')
hold on
stem([0:N/2-1],a_exact(1:N/2),'or')
hold off
axis([0,N/2,-0.15,1.3])
xlabel('n')
ylabel('a_n')
legend('Computed a_n','Exact a_n')
title('Plot of computed and exact values of a_n')

% figure(2)
% stem([0:N-1],real(z))
% title('Subplot 1: real(z)')
% 
% figure(3)
% stem([0:N-1],imag(z));
% title('Subplot 2: imag(z)')
% 
% figure(4)
% plot(x,y);
% title('Subplot 4: g(x)')

%% assignment 3 (a)
clear all;
clc;

h = @(x) 1.*((0<x) & (x < pi)) + (-1).*((pi < x) & (x < 2*pi));

N = 2^8;
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
x_5 = part_sum_sincos(z,x2,5);
x_10 = part_sum_sincos(z,x2,10);
x_50 = part_sum_sincos(z,x2,50);

figure(1);
stem([0:N-1],real(z));
hold on;
stem([0:N-1],imag(z));

figure(2);
plot(x2,h(x2), '-b');
hold on;
plot(x2,x_5,'-r');
hold on;
plot(x2,x_10, '-g');
hold on;
plot(x2,x_50, '-c');
legend('h(x)','Partial sum M=5','Partial sum M=10', 'Partial sum M=50');
xlabel('x');
ylabel('y');
title('h(x) together with fourier partial sum of h(x), N=2^8');


%% assignment 3 (b)
clear all;
clc;
h = @(x) 1.*((0<x) & (x < pi)) + (-1).*((pi < x) & (x < 2*pi));

N = 2^16;
x = zeros(1,N);
% (0,2pi)
for i=0:(N-1)
    x(i+1) = (2*pi*i/N);
end
y = h(x);
z=mydft(y);
%z = fft(y)/N;
x2 = 0:0.001:2*pi;

M = [1,5,10,25,50,100,200,500,1000,2000];
%M = [1000,2000,5000,10000,20000,30000];
E = zeros(1,length(M));

for i=1:length(M)
   E(i) = part_sum_max_error(z,x2,h,M(i));
end

figure(1);
plot(M,E,'-r');
legend('Error max(h(x)-s_M(x))');
xlabel('M');
ylabel('Error');
title('Error as a function of number of terms in partial sums');

%% assignment 4 (a)
clear all;
clc;

runsPerN = 5;
N = [2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10, 2^11, 2^12];
AvgDFTTimes = zeros(1,length(N));
AvgFFTTimes = zeros(1,length(N));



for i = 1:length(N)
    
    for j = 1:runsPerN
        vector = rand(1,N(i));
        result = zeros(1,N(i));
        
        tic;
        result = mydft(vector);
        AvgDFTTimes(i) = AvgDFTTimes(i) + toc;
        
        tic;
        result = fft(vector)/N(i);
        AvgFFTTimes(i) = AvgFFTTimes(i) + toc;
        
    end
    AvgDFTTimes(i) = AvgDFTTimes(i)/runsPerN;
    AvgFFTTimes(i) = AvgFFTTimes(i)/runsPerN;
    disp(N(i));
end

figure(1);
plot(N,AvgDFTTimes,'-r');
hold on;
plot(N,AvgFFTTimes,'-b');
legend('Average DFT Time', 'Average FFT Time');
xlabel('Vector length (N)');
ylabel('Time (s)');
title('Running time for DFT versus FFT');

TimeN12 = AvgDFTTimes(end) - AvgFFTTimes(end);

%avgDFTTime / avgFFTTime


%% assignment 4 (b)
clear all;
clc;

fid = fopen('filtre.data.txt','r');
Y = fscanf(fid,'%f',[1,inf]);
fclose(fid);

Z = fft(Y); 

cutoff = 0.5;

toRemove = find(abs(Z) < cutoff);

Z(toRemove) = 0;

Y_new = ifft(Z);

x2 = 0:(length(Y)-1);
x2 = x2*2*pi/(length(Y)-1);
figure(1);
plot(x2,Y,'-b');
legend('Data');
xlabel('x');
ylabel('y');
title('Data');
axis([0,2*pi,-1.1,1.1]);
figure(2);
plot(x2,Y_new, '-r');
legend('Filtered Data');
xlabel('x');
ylabel('y');
title('Filtered data');
axis([0,2*pi,-1.1,1.1]);


%% assignment 4 (b)
clear all;
clc;

%load('handel');
load('splat');
%load('gong');
sound(y,Fs);

N = length(y);
Z = fft(y);
W = Z;
M = max(abs(Z));
cutoff = 2;

toRemove = find(abs(W) < cutoff);
W(toRemove) = 0;

ZS = sparse(Z);
WS = sparse(W);

before = whos('ZS');
after = whos('WS');
comprRatio = before.bytes/after.bytes

pause(10);
disp('Play compressed signal');
w = real(ifft(full(WS)));
sound(w,Fs);

