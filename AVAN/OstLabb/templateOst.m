%%%%%%%%%%%% INFORMATION %%%%%%%%%%%%%%%%


%%%%%%%%%%%%% GLOBAL DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%
%%%% GENERAL SETTINGS %%%%
dt = 0.001;
drawUpdate = 2; %update every drawUpdate:th frame
energyUpdate = 1; %update every energyUpdate:th frame

TIMES = 0:dt:4;

DIM = 2; %decides how to draw, setup ground and in which direction gravity acts
G = -5; %gravity constant

%%%% PARTICLES SETTINGS %%%%
particleColumns = 7; %number of columns x-axis for particle system
particleRows = 4; %number of rows y-axis for particle system
particleLayers = 3; %number of layers z-axis for particle system
if DIM == 2
   particleLayers = 1;
end
particleNumbers = particleColumns*particleRows*particleLayers; %number of particles
particleScale = 1; %particle-scale
particleMasses = ones(particleNumbers,1)*1;

%%%% springs SETTINGS %%%%
springNumbers = get_number_of_springs(particleColumns,particleRows,particleLayers);
springKS = 200; %coefficient for spring-force
springKD = 0; %coefficient of spring-damping

%%%% ENERGY SETTINGS %%%%
energyCalcN = floor(length(TIMES)/energyUpdate)+1;
energyTimes = zeros(energyCalcN, 1);

particleEnergyKinetic = zeros(energyCalcN, 1);
particleEnergyPotential = zeros(energyCalcN, 1);
particleEnergySprings = zeros(energyCalcN, 1);
particleEnergyTotal = zeros(energyCalcN, 1);

%%%% GROUND SETTINGS %%%%
groundNX = 0;
groundNY = 0;
if DIM == 2
    groundNX = 30;
    groundNY = 20;
elseif DIM == 3
    groundNX = 30; %number of balls in x-dir
    groundNY = 15; %number of balls in y-dir 
end
groundScale = 0.5; %ground-scale
groundStartPoint = [0,0,0]; %ground-starting-point
groundBallRadius=0.5;  % Radius for ground balls

%%%%%%%%%%% GENERATE GROUND %%%%%%%%%%%%%%
groundBallPositions = generate_ground(groundNX,groundNY,groundScale,groundStartPoint,groundBallRadius,DIM);
%%% SETUP BALLS GRAPHICAL OBJECTS %%%

%%%%%%%%%% GENERATE CUBOID POSITIONS AND springs %%%%%%%
[particlePositions, springs] = setup_cuboid(particleColumns,particleRows,particleLayers,particleScale);
%%%% SETUP CUBOID lines %%%%
lines = setup_lines_cross(particlePositions, springs, DIM);

%%%%%%%%% PARTICLE PERTRUBATIONS %%%%%%%%
%%%% INITIAL POSITION %%%%
particlePositionsOld = particlePositions;
if DIM == 2
    particlePositions(:,:) = particlePositions(:,:) + [2,2,0];
elseif DIM == 3
    particlePositions(:,:) = particlePositions(:,:) + [0,0,2];
end
%%%% INITIAL VELOCITIES %%%%
particleVelocities = zeros(particleNumbers,3);
particleVelocitiesOld = particleVelocities;
if DIM == 2
    particleVelocities(:,:) = particleVelocities(:,:) + [0,0,0];
elseif DIM == 3
    particleVelocities(:,:) = particleVelocities(:,:) + [2,0,0];
end
%%%% INITIAL FORCES %%%%
particleForces = zeros(particleNumbers, 3);

%%%%%%%%%%%%%%% GRAPHICS SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
%figure(1);
if DIM == 2
    axis equal; 
elseif DIM == 3
    daspect([1,1,1]);
    axis([...
        -5*groundBallRadius*groundScale+groundStartPoint(1),...
        groundNX*groundScale+5*groundBallRadius*groundScale+groundStartPoint(1),...
        -5*groundBallRadius*groundScale+groundStartPoint(2),...
        groundNY*groundScale+5*groundBallRadius*groundScale+groundStartPoint(2),...
        groundStartPoint(3),...
        groundStartPoint(3)+10]...
        );
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
end

%%%% INITIAL ENERGY UPDATE %%%%
particleEnergyKinetic(1) = energy_kinetic(particleVelocities, particleMasses);
particleEnergySprings(1) = energy_springs(particlePositions,springs,springKS);
particleEnergyPotential(1) = energy_gravity(particlePositions,particleMasses,G,DIM);
particleEnergyTotal(1) = particleEnergyKinetic(1) + particleEnergySprings(1) + particleEnergyPotential(1);

%%%% INITIAL UPDATE %%%%
[particleVelocities,particleForces] = init_update(dt,particlePositions,particleVelocities,G,particleMasses,springKS,springKD,springs,DIM);

%profile on;

%%%% MAIN LOOP %%%%
atEnergies = 2;
updateNumber = 1;
for t = TIMES
    %half euler step forward to calculate forces in sync
    particleVelocitiesHalfForward = particleVelocities + 0.5*dt*particleForces./particleMasses;
    
    particleForces = zeros(particleNumbers,3);
    %UPDATE GRAVITY
    particleForces = calc_gravity_forces(particleForces,particleMasses,G,DIM);
    %UPDATE SPRING FORCES, use half euler forward
    particleForces = calc_spring_forces(particlePositions,particleVelocitiesHalfForward,particleForces,springKS,springKD,springs);
    %UPDATE POSITION AND VELOCITY
    [particlePositions, particleVelocities] = update_rv(dt,particlePositions,particleVelocities,particleForces,particleMasses);
    
    %UPDATE COLLISION
    [particlePositions, particleVelocities] = update_coll(particlePositions,particleVelocities,groundBallPositions,groundBallRadius);
    
    % take a half euler step forward to calculate energies in sync
    particleVelocitiesHalfForward = particleVelocities + 0.5*dt*particleForces./particleMasses;
    
    if mod(updateNumber,drawUpdate) == 0
        lines = draw_lines(particlePositions, lines);
        drawnow;
    end
    
    %%%% PLOT AND CALC ENERGIES %%%%
    
    if mod(updateNumber,energyUpdate) == 0
        energyTimes(atEnergies) = t;
        % PARTICLE ENERGIES
        particleEnergyKinetic(atEnergies) = energy_kinetic(particleVelocitiesHalfForward, particleMasses);
        particleEnergySprings(atEnergies) = energy_springs(particlePositions,springs,springKS);
        particleEnergyPotential(atEnergies) = energy_gravity(particlePositions,particleMasses,G,DIM);
        particleEnergyTotal(atEnergies) = particleEnergyKinetic(atEnergies)...
            + particleEnergySprings(atEnergies)...
            + particleEnergyPotential(atEnergies);
        
        atEnergies = atEnergies + 1;
    end
    
    updateNumber = updateNumber + 1;
end

%profile viewer;

figure(2);
plot(energyTimes, particleEnergyKinetic, '-r');
hold on;
plot(energyTimes, particleEnergySprings, '-b');
hold on;
plot(energyTimes, particleEnergyPotential, '-g');
hold on;
plot(energyTimes, particleEnergyTotal, '-m');
legend('Kinetic','springs','Potential','Total');

maxEnergyDifference = max_energy_diff(particleEnergyTotal)






