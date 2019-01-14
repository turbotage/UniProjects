%%%%%%%%%%%% INFORMATION %%%%%%%%%%%%%%%%
% This is an example of excersize 2 without gravity or ground
% Whatch spring effects by applying changes in the perturbations section

%%%%%%%%%%%%% GLOBAL DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%
dt = 0.005;
drawUpdate = dt*10;
plotUpdate = dt*10;
energyUpdate = dt*10;

TIMES = 0:dt:1000;

G = -1; %gravity constant

Ek = zeros(length(TIMES), 1);
Ep = zeros(length(TIMES), 1);
Ef = zeros(length(TIMES), 1);
E = zeros(length(TIMES), 1);

SP = [1,5,1]; % Start-position (centrum)
Ra=0.7;  % Radie

NR = 4;
NC = 8;
NP = NR*NC;
NS = 4*NR*NC - 3*NR - 3*NC + 2;
Ks = 70;
Kd = 2;
Kf = 0.5;
M = ones(NP,1);
dim = 3;

%%%%%%%%%%%%%%% GRAPHICS SETTINGS %%%%%%%%%%%%%%%%%%%%%%%% 
daspect([1,1,1]);   % Skalar axlarna lika
% Axlarnas intervall
if dim == 3
    axis([0,10,0,10,0,3]);
elseif dim == 2
    %axis([0,40,-20,20]);  % stairs
    axis([0, 20, 0, 20]);
end
hold on;             % Laas fast dessa instaellningar


% number of rows, number of cols, scale, initial translation
%%%%%%%%%% GENERATE QUAD POSITIONS %%%%%%%
PR_n = generate_quad(NR,NC,1,SP);
PR_n(:,2) = PR_n(:,2) + ones(NP,1);
%%%% INITIAL VELOCITY %%%%
PV_n = zeros(NP,3);
PV_n(:,1) = ones(NP,1)*3;
%%%% INITIAL FORCE %%%%%%
PF_n = zeros(NP,3);




%%%%%%%% SETUP SPRINGS FROM POSITIONS %%%%%%%%%%%
springs = setup_springs(PR_n, NR, NC, Ks, Kd);
%lines = setup_lines(PR_n, springs, NR, NC, dim);
lines = setup_lines_cross(PR_n, springs, dim);




%%%%%%%%%%% GENERATE GROUND %%%%%%%%%%%%%%
[BR_n, r_n] = generate_ground_stairs(8,4,10,1, [0,5,1]);
%[BR_n, r_n] = generate_ground(40,0.5,SP);
ground_balls = setup_balls(BR_n, r_n, 2);




%%%%%%%%%% PERTURBATIONS %%%%%%%%%%%%%
%for i=1:NC %perturb upper row
%    PV_n(i,:) = PV_n(i,:) + [3, 0, 0];
%end
%perI = 1; %change to change what particle that gets perturbed
%PR_n(perI,:) = PR_n(perI,:) + [0.1, 0.1, 0];
%PV_n(perI,:) = PV_n(perI,:) + [0, 0, 0];
%PF_n(perI,:) = PF_n(perI,:) + [1, 1, 0];




PV_n = init_update(dt,PR_n,PV_n,PF_n,springs,M);
%%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%
for t = 0:dt:1000
    PF_n = zeros(NP,3);
    %UPDATE GRAVITY
    PF_n = update_g(PF_n, M, G, dim);
    %UPDATE COLLISION
    [PR_n, PV_n] = update_coll(PR_n, PV_n, ground_balls, BR_n, Kf);
    %UPDATE POSITION AND VELOCITY
    [PR_n, PV_n] = update_rv(dt, PR_n, PV_n, PF_n, springs, M);
    
    if mod(t,drawUpdate) == 0
        lines = draw_lines(PR_n, lines);
        drawnow;
    end
end