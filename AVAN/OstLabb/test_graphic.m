xa=2;    % Start-position (centrum)
ya=5;
za=1;
Da=0;    % Vinkeln
Ra=0.7;  % Radie
La=0.4;  % Kryssets storlek

NR = 2;
NC = 2;
NP = NR*NC;
NS = 4*NR*NC - 3*NR - 3*NC + 2;
Ks = 2;
Kd = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instaellningar for grafik-fonstret
daspect([1,1,1]);   % Skalar axlarna lika
axis([0,8,0,6]);    % Axlarnas intervall
hold on;             % Laas fast dessa instaellningar


% SETUP POSITIONS
%R_n = [1,3,1; 2,3,1; 3,3,1; 1,2,1; 2,2,1; 3,2,1; 1,1,1; 2,1,1; 3,1,1];

% number of rows, number of cols, scale, initial translation
R_n = generate_quad(NR,NC,0.4, [0,1,0]);
V_n = [];
F_n = [];

% SETUP SPRINGS FROM POSITIONS
springs = setup_springs(R_n, NR, NC, Ks, Kd);





