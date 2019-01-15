%%%%%%%%%%%% INFORMATION %%%%%%%%%%%%%%%%
% This is an example of excersize 2 without gravity or ground
% Whatch spring effects by applying changes in the perturbations section

%%%%%%%%%%%%% GLOBAL DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%
dt = 0.005;
drawUpdate = 10;
energyUpdate = 10;

TIMES = 0:dt:2;

G = -5; %gravity constant

energyN = (length(TIMES)-1)/energyUpdate;
Etimes = zeros(energyN, 1);
Ek = zeros(energyN, 1);
Ep = zeros(energyN, 1);
Ef = zeros(energyN, 1);
E = zeros(energyN, 1);

SP = [0,5,1]; % Start-position (centrum)
Ra=0.7;  % Radie

NR = 4;
NC = 8;
NL = 1;
NP = NR*NC;
NS = 4*NR*NC - 3*NR - 3*NC + 2;
Ks = 70;
Kd = 2;
Kf = 0.5;
M = ones(NP,1);
dim = 2;

%%%%%%%%%%%%%%% GRAPHICS SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
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
%%%%%%%%%% GENERATE QUAD POSITIONS AND SPRINGS %%%%%%%
PR_n = generate_quad(NR,NC,1,SP + [0,1,0]);
%%%% INITIAL VELOCITY %%%%
PV_n = zeros(NP,3);
PV_n(:,1) = ones(NP,1)*8;
%%%% INITIAL FORCE %%%%%%
PF_n = zeros(NP,3);

%%%% SPRINGS %%%%%
springs = setup_springs(PR_n, NR, NC, Ks, Kd);

%%%%%%%%%%% GENERATE GROUND %%%%%%%%%%%%%%
%[BR_n, r_n] = generate_ground_stairs(8,4,10,1, [0,5,1]);
[BR_n, r_n] = generate_ground(40,0.5,SP);

%%%%%%%%%%%%%% SETUP GRAPHICAL OBJECTS %%%%%%%%%%%%%
%LINES
%lines = setup_lines(PR_n, springs, NR, NC, dim);
lines = setup_lines_cross(PR_n, springs, dim);
%BALLS
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
atEnergies = 1;
%%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%
for t = TIMES
    PF_n = zeros(NP,3);
    %UPDATE GRAVITY
    PF_n = update_g(PF_n, M, G, dim);
    %UPDATE COLLISION
    [PR_n, PV_n] = update_coll(PR_n, PV_n, ground_balls, BR_n, Kf);
    %UPDATE POSITION AND VELOCITY
    [PR_n, PV_n] = update_rv(dt, PR_n, PV_n, PF_n, springs, M);
    
    if mod(t,dt*drawUpdate) == 0
        lines = draw_lines(PR_n, lines);
        drawnow;
    end
    
    if mod(t,dt*energyUpdate) == 0
        Etimes(atEnergies) = t;
        Ek(atEnergies) = energy_kinetic(PV_n, M);
        Ef(atEnergies) = energy_springs(PR_n,springs);
        Ep(atEnergies) = energy_gravity(PR_n,M,G,dim);
        E(atEnergies) = Ek(atEnergies) + Ef(atEnergies) + Ep(atEnergies);
        atEnergies = atEnergies + 1;
    end
end

figure(2);
plot(Etimes, Ek, '-r');
hold on;
plot(Etimes, Ef, '-b');
hold on;
plot(Etimes, Ep, '-g');
hold on;
plot(Etimes, E, '-m');

