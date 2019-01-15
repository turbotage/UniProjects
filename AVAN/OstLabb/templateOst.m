%%%%%%%%%%%% INFORMATION %%%%%%%%%%%%%%%%


%%%%%%%%%%%%% GLOBAL DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%
%%%% GENERAL SETTINGS %%%%
dt = 0.004;
drawUpdate = 10; %update every drawUpdate:th frame
energyUpdate = 10; %update every energyUpdate:th frame

TIMES = 0:dt:8;

DIM = 3; %decides how to draw, getup ground and in which direction gravity acts
G = -5; %gravity constant

%%%% PARTICLES SETTINGS %%%%
NC = 3; %number of columns x-axis for particle system
NR = 3; %number of rows y-axis for particle system
NL = 4; %number of layers z-axis for particle system
if DIM == 2
   NL = 1; 
end
NP = NR*NC*NL; %number of particles
NS = NL*NC*(NR-1)+NR*NC*(NL-1)+NL*NR*(NC-1)+... %straights
     2*NL*(NR-1)*(NC-1)+2*NR*(NC-1)*(NL-1)+2*NC*(NL-1)*(NR-1)+... %2cross
     4*(NL-1)*(NR-1)*(NC-1); %3cross
KS = 200; %coefficient for spring-force
KD = 0; %coefficient of spring-damping
PSP = [0,0,0]; %particle-start-point
PS = 1; %particle-scale
PM = ones(NP,1)*0.5;

EN = (length(TIMES)-1)/energyUpdate;
ETIMES = zeros(EN, 1);
EKP = zeros(EN, 1);
EPP = zeros(EN, 1);
EFP = zeros(EN, 1);
ETP = zeros(EN, 1);


%%%% GROUND SETTINGS %%%%
NX = 20; %number of balls in x-dir
NY = 20; %number of balls in y-dir
GS = 0.25; %ground-scale
GSP = [0,0,0]; %ground-starting-point
GBR=0.7;  % Radius for ground balls

%%%%%%%%%%%%%%% GRAPHICS SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
if DIM == 2
    axis equal; 
elseif DIM == 3
    axis([0,NX*GBR,0,NX*GBR,0,NX*GBR]);
end

%%%%%%%%%%% GENERATE GROUND %%%%%%%%%%%%%%
BR_n = generate_ground(NX,NY,GS,GSP,GBR,DIM);
%%% SETUP BALLS GRAPHICAL OBJECTS %%%

%%%%%%%%%% GENERATE CUBOID POSITIONS AND SPRINGS %%%%%%%
[PR_n, springs] = setup_cuboid(NC,NR,NL,PS);
%%%% SETUP CUBOID LINES %%%%
lines = setup_lines_cross(PR_n, springs, DIM);

%%%%%%%%% PARTICLE PERTRUBATIONS %%%%%%%%
%%%% INITIAL POSITION %%%%
if DIM == 2
    PR_n(:,:) = PR_n(:,:) + [2,2,0];
elseif DIM == 3
    PR_n(:,3) = PR_n(:,3) + ones(NP,1)*3;
end
%%%% INITIAL VELOCITIES %%%%
PV_n = zeros(NP,3);
if DIM == 2
    PV_n(:,:) = PV_n(:,:) + [3,0,0];
elseif DIM == 3
    PV_n(:,:) = PV_n(:,:) + [2,1,0];
end

%%%% INITIAL UPDATE %%%%
[PV_n,PF_n] = init_update(dt,PR_n,PV_n,G,PM,KS,KD,springs,DIM);

%%%% MAIN LOOP %%%%
atEnergies = 1;
for t = TIMES
    %half euler step forward
    V_hf = PV_n + 0.5*dt*PF_n./PM;
    
    %UPDATE COLLISION
    [PR_n, PV_n] = update_coll(PR_n,PV_n,BR_n,GBR);
    
    %%%% BIG CALC
    PF_n = zeros(NP,3);
    %UPDATE GRAVITY
    PF_n = calc_gravity_forces(PF_n,PM,G,DIM);
    %UPDATE SPRING FORCES, use half euler forward
    PF_n = calc_spring_forces(PR_n,V_hf,PF_n,KS,KD,springs);
    %UPDATE POSITION AND VELOCITY
    [PR_n, PV_n] = update_rv(dt,PR_n,PV_n,PF_n,PM);
    
    if mod(t,dt*drawUpdate) == 0
        lines = draw_lines(PR_n, lines);
        drawnow;
    end
    
    %%%% PLOT AND CALC ENERGIES %%%%
    
    if mod(t,dt*energyUpdate) == 0
        ETIMES(atEnergies) = t;
        % PARTICLE ENERGIES
        EKP(atEnergies) = energy_kinetic(PV_n, PM);
        EFP(atEnergies) = energy_springs(PR_n,springs,KS);
        EPP(atEnergies) = energy_gravity(PR_n,PM,G,DIM);
        ETP(atEnergies) = EKP(atEnergies) + EFP(atEnergies) + EPP(atEnergies);
        
        atEnergies = atEnergies + 1;
    end
end

figure(2);
plot(ETIMES, EKP, '-r');
hold on;
plot(ETIMES, EFP, '-b');
hold on;
plot(ETIMES, EPP, '-g');
hold on;
plot(ETIMES, ETP, '-m');
legend('Kinetic','Springs','Potential','Total');






