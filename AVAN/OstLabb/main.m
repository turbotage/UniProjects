%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global definitians
dt = 0.001;
drawGeomN = dt*10;
calcEnergies = dt*10;

TIMES = 0:dt:1000;

xa=1;    % Start-position (centrum)
ya=1;
za=1;
Da=0;    % Vinkeln
Ra=0.7;  % Radie
La=0.4;  % Kryssets storlek
SP = [xa,ya,za];

NR = 1;
NC = 2;
NP = NR*NC;
NS = 4*NR*NC - 3*NR - 3*NC + 2;
Ks = 100;
Kd = 0.5;
M = ones(NP,1);
dim = 2;

energyN = length(TIMES)/calcEnergies;
Ek = zeros(energyN, 1);
Ep = zeros(energyN, 1);
Ef = zeros(energyN, 1);
E = zeros(energyN, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for 
figure(1);
daspect([1,1,1]);   % Skalar axlarna lika
% Axlarnas intervall
if dim == 3
    axis([0,10,0,10,0,3]);
elseif dim == 2
    axis([0,10,0,10]);
end
hold on;             % Laas fast dessa instaellningar


% SETUP POSITIONS
%R_n = [1,3,1; 2,3,1; 3,3,1; 1,2,1; 2,2,1; 3,2,1; 1,1,1; 2,1,1; 3,1,1];

% number of rows, number of cols, scale, initial translation
R_n = generate_quad(NR,NC,1, SP);
V_n = zeros(NP,3);
F_n = zeros(NP,3);

% SETUP SPRINGS FROM POSITIONS
springs = setup_springs(R_n, NR, NC, Ks, Kd);
%lines = setup_lines(R_n, springs, NR, NC, dim);
lines = setup_lines_cross(R_n, springs, dim);
%balls = setup_balls(lines);

% offset one particle 
R_n(1,:) = R_n(1,:) + [0.5, 0.5, 0];
%F_n(1,:) = F_n(1,:) + [1, 1, 0]; 


V_n = init_update(dt,R_n,V_n,F_n,springs,M);

%MAIN LOOP
for t = TIMES
    F_n = zeros(NP,3);
    
    [R_n, V_n] = update_rv(dt,R_n,V_n,F_n,springs,M);
    if mod(t,drawEveryN) == 0 % draw geometry
        figure(1);
        lines = draw_lines(R_n, lines);
        drawnow;
    end
    %if mod(t,drawPlotN) == 0
    %    figure(2);
    %end
end