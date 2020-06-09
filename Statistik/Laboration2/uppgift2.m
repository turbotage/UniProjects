% load and check assumptions
clc;
clear all;

IN = load('kroppstemperatur.txt');

T = IN(:,1);
X = IN(:,3);

nmen = 0;
nwomen = 0;

for i=1:length(T)
    if IN(i,2) == 1
        nmen = nmen + 1;
    else
        nwomen = nwomen + 1;
    end
end

T_m = T(1:nmen);
T_k = T(nmen+1:length(T));



% assumption 1
[h_T,p_T] = lillietest(T)
[h_X,p_T] = lillietest(X)

% assumption 2
[h_T_m,p_T_m] = lillietest(T_m)
[h_T_k,p_T_k] = lillietest(T_k)

[h_v,p_v] = vartest2(T_m,T_k)


%%
%Problem 1
[h,p] = ttest(T,98.6)

%%
%Problem 2

[h,p] = ttest2(T_m,T_k)
	
%%
%Problem 3
T_reg = [ones(size(T)),T];
[B,BINT,R,RINT,STATS] = regress(X,T_reg)

f = @(x) B(2).*x + B(1)

figure(1);
plot(T,X,'.');
ylabel('Hjärtfrekvens [BPM]');
xlabel('Kroppstemperatur [F]');
title('Spridningsdiagram av datat');

figure(2);
plot(T,X,'.');
hold on;
fplot(f);
axis([96.2 100.9 56 90]);
ylabel('Hjärtfrekvens [BPM]');
xlabel('Kroppstemperatur [F]');
title('Spridningsdiagram och regressionslinje');


