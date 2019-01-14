xa=1;    % Start-position (centrum)
ya=1;
za=1;
Da=0;    % Vinkeln
Ra=0.7;  % Radie
La=0.4;  % Kryssets storlek
SP = [xa,ya,za];

NR = 4;
NC = 4;
NP = NR*NC;
NS = 4*NR*NC - 3*NR - 3*NC + 2;
Ks = 100;
Kd = 0.5;
M = ones(NP,1);
A = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instaellningar for grafik-fonstret
daspect([1,1,1]);   % Skalar axlarna lika
axis([0,8,0,6]);    % Axlarnas intervall
hold on;             % Laas fast dessa instaellningar


% SETUP POSITIONS
%R_n = [1,3,1; 2,3,1; 3,3,1; 1,2,1; 2,2,1; 3,2,1; 1,1,1; 2,1,1; 3,1,1];

% number of rows, number of cols, scale, initial translation
R_n = generate_quad(NR,NC,1, SP);
V_n = zeros(NP,3);
F_n = zeros(NP,3);

% SETUP SPRINGS FROM POSITIONS
springs = setup_springs(R_n, NR, NC, Ks, Kd);
lines = setup_lines_cross(R_n, springs, NR, NC);
%lines = setup_lines_cross(R_n, springs, NR, NC);
%balls = setup_balls(lines);

% offset one particle 
R_n(1,:) = R_n(1,:) + [0.1, 0.1, 0.1];
%F_n(1,:) = F_n(1,:) + [1, 1, 0]; 


dt = 0.001;
drawEveryN = dt*10;

V_n = init_update(dt,R_n,V_n,F_n,springs,M);

for t = 0:dt:1000
    [R_n, V_n] = update_RandV(dt,R_n,V_n,springs,M);
    if mod(t,drawEveryN) == 0
        lines = draw(R_n, lines);
        drawnow;
    end
end






