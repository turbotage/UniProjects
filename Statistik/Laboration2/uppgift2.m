%%
%uppgift 1
T = load('kroppstemperatur.txt');

[h] = lillietest(T(:,1)) % det går inte att förkasta att temperaturen är normalfördelad

h = ttest(T(:,1),98.6) % ger h=1 dvs h_0 kan förkastas

%%
%uppgift 2

nmen = 0;
nwomen = 0;

for i=1:length(T(:,1))
    if T(i,2) == 1
        nmen = nmen + 1;
    else
        nwomen = nwomen + 1;
    end
end

menTemperatures = T(1:nmen,1);
womenTemperatures = T(nmen+1:length(T(:,1)),1);

h_men = lillietest(menTemperatures)
h_women = lillietest(womenTemperatures)

h = ttest2(menTemperatures, womenTemperatures)
	
%%
%uppgift3

x = T(:,1);
y = T(:,3);


nmen = 0;
nwomen = 0;

for i=1:length(T(:,1))
    if T(i,2) == 1
        nmen = nmen + 1;
    else
        nwomen = nwomen + 1;
    end
end

menTemperatures = T(1:nmen,1);
menHeartF = T(1:nmen,3);
womenTemperatures = T(nmen+1:length(T(:,1)),1);
womenHeartF = T(nmen+1:length(T(:,1)),3);


X = [ones(size(x)),x];
[B,BINT,R,RINT,STATS] = regress(y,X)

f = @(x) B(2)*x + B(1)

figure(1);
plot(T(:,1),T(:,3),'.');
hold on;
fplot(f);
axis auto xy;


