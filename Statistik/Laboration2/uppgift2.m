T = load('kroppstemperatur.txt');

at = mean(T(:,1)); %average temperature

atdev = std(T); % standard deviation of average temperature

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

ment = T(1:nmen,1);
woment = T(nmen+1:length(T(:,1)),1);

menavgt = mean(ment);
womenavgt = mean(woment);

menstdt = std(ment);
womenstdt = std(woment);

%uppgift3

x = T(:,1);
y = T(:,3);
X = [ones(size(x)),x];
[B,BINT,R,RINT,STATS] = regress(y,X)

f = @(x) B(2)*x + B(1)

figure(1)
plot(x,y,'.');
hold on;
fplot(f);
axis([95 102 50 100]);


