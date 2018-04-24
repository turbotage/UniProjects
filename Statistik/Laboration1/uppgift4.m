
% a)

x = poissrnd(5, 30, 1000);

means = size(30,1);

for i = 1:1000
   means(i) = mean(x(:,i));
end

% b)

hist(means,1000)

% c)

mean(means)
var(means)