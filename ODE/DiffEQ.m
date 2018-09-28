
%% 
%assignment1
clear all;
clc; 
load('randwebs.mat');
load('randinits.mat');

intmat = 5;

options = odeset('NonNegative',1);

A = randwebscell{intmat};
tend = 1000;
tspan = [0,tend];
y0 = randinitscell{intmat};

f = @(t,x) x.*A*x;


[t,x] = ode45(f,tspan,y0,options);
plot(t,x)
hold on
xlabel('Time, t')
hold on
ylabel('Quantity of Specis, UNIT???')
%% 
%assignment2
clear all;
clc; 
load('randwebs.mat');
load('randinits.mat');

intmat = 1;

options = odeset('NonNegative',1);

A = randwebscell{intmat};
tend = 2000;
tspan = [0,tend];
y0 = randinitscell{intmat};

f = @(t,x) x.*A*x

[t,x] = ode45(f,tspan,y0,options);
plot(t,x)
hold on
xlabel('Time, t')
hold on
ylabel('Quantity of Specis, UNIT???')

aliveThreshold = exp(-5)
xalive = zeros(1,length(x(1,:)));

for i = (length(x(:,1))-2):length(x(:,1))
    xalive = xalive + ( x(i,:) > aliveThreshold );
end
xalive = (xalive > 0);

display(xalive);
sum(xalive)

%% 
%assignment3and4 and optional plot for assignment 5
clear all;
clc; 
load('randwebs.mat');
load('randinits.mat');

options = odeset('NonNegative',1);

tend = 10000;

interactionset_dead = zeros(1000,1);
species_alive = zeros(11,1);

hasPlotted = 0;

for k = 1:1000

    intmat = k;

    A = randwebscell{intmat};
    tspan = [0,tend];
    y0 = randinitscell{intmat};

    f = @(t,x) x.*A*x;

    [t,x] = ode45(f,tspan,y0,options);
    
    %change for different threshhold values
    aliveThreshold = exp(-5)
    xalive = zeros(1,length(x(1,:)));

    for i = (length(x(:,1))-2):length(x(:,1))
        xalive = xalive + ( x(i,:) > aliveThreshold );
    end
    xalive = (xalive > 0);
    
    alivesum = sum(xalive);
    interactionset_dead(k) = alivesum;
    species_alive(alivesum+1) = species_alive(alivesum+1) + 1;
    
    %uncommnent to se plot of all species alive, assignment 5
    %{
    if (s==10) && (hasPlotted == 0)
        plot(t,x)
        hold on
        xlabel('Time, t')
        hold on
        ylabel('Quantity of Specis')
        hasPlotted = 1;
    end
    %}
    
end

display(interactionset_dead); %how many species that die in each set
display(species_alive) % how many times 0,1..10 species stays alive 

%% 
%assignment9
clear all;
clc; 
load('randwebs.mat');
load('randinits.mat');

%format LONGE;

options = odeset('NonNegative',1);

%changes for how long we solves the equations
tend = 10000;
tspan = [0,tend];

%change to check which microbes are oscilating in corresponding
%initialconditions and to plot the solution
initialConditionsNum = 1000;

for k = 1:1000

    intmat = k;
    
    A = randwebscell{intmat};
    y0 = randinitscell{intmat};

    f = @(t,x) x.*A*x;

    [t,x] = ode45(f,tspan,y0,options);
    
    smallThreshold = exp(-15);
    start = length(x(:,1)) - floor(length(x(:,1))*0.15);
    checkDerivativeSignArray = x(start:end,:);
    signChanges = zeros(1,10);
    Y = diff(checkDerivativeSignArray);
    for i=1:(length(Y(:,1))-1)
       temp = ( (Y(i,:).*Y(i+1,:)) < 0);
       temp2 = (Y(i,:) > smallThreshold); % don't consider signchange for small values since we can have +-0
       signChanges = signChanges + (temp.*temp2);
    end
    
    if k == initialConditionsNum
        plot(t,x)
        hold on
        xlabel('Time, t')
        hold on
        ylabel('Quantity of Specis')
        disp(signChanges); %shows how many times the derivative makes a signchange
        disp(signChanges > 0); %shows how many microbes that oscilates in the end of the solution
    end
    
end

