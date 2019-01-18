%%%%%%%%%%%% INFORMATION %%%%%%%%%%%%%%%%


%%%%%%%%%%%%% GLOBAL DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%
%%%% GENERAL SETTINGS %%%%
dt = 0.001;
drawUpdate = 5; %update every drawUpdate:th frame
energyUpdate = 20; %update every energyUpdate:th frame

TIMES = 0:dt:1;

DIM = 3; %decides how to draw, setup ground and in which direction gravity acts
G = -1; %gravity constant

%%%% PARTICLES SETTINGS %%%%
particleColumns = 12; %number of columns x-axis for particle system
particleRows = 12; %number of rows y-axis for particle system
particleLayers = 1; %number of layers z-axis for particle system
particleNumbers = particleColumns*particleRows*particleLayers; %number of particles
particleScale = 0.4; %particle-scale
particleMasses = ones(particleNumbers,1)*0.1;

%%%% springs SETTINGS %%%%
springNumbers = get_number_of_springs(particleColumns,particleRows,particleLayers);
springKS = 10000; %coefficient for spring-force
springKD = 20; %coefficient of spring-damping

%%%% ENERGY SETTINGS %%%%
energyCalcN = floor(length(TIMES)/energyUpdate)+1;
energyTimes = zeros(energyCalcN, 1);

particleEnergyKinetic = zeros(energyCalcN, 1);
particleEnergyPotential = zeros(energyCalcN, 1);
particleEnergySprings = zeros(energyCalcN, 1);
particleEnergyTotal = zeros(energyCalcN, 1);

ballEnergyKinetic = zeros(energyCalcN,1);
ballEnergyPotential = zeros(energyCalcN,1);
ballEnergyTotal = zeros(energyCalcN,1);

systemEnergyTotal = zeros(energyCalcN,1);

%%%%%%%%%% GENERATE BALL AND BALL SETTINGS %%%%%%%%%%
ballN = 1;
ballMasses = ones(ballN,1)*20;
ballPositions = [2,2,6];
ballVelocities = [0,0,0];
ballForces = [0,0,0];
ballRadiuses = [0.8];
[bx,by,bz] = sphere(15);
balls = [];
for i=1:length(ballPositions(:,1))
   balls(i) = surf(bx*ballRadiuses(i) + ballPositions(1,1), by*ballRadiuses(i) + ballPositions(1,2), bz*ballRadiuses(i) + ballPositions(1,3), 'FaceColor', 'r'); 
end

%%%%%%%%%% GENERATE CUBOID POSITIONS AND springs %%%%%%%
[particlePositions, springs] = setup_cuboid(particleColumns,particleRows,particleLayers,particleScale);
particleInnerIndices = get_inner_indices(particleColumns,particleRows,particleLayers);
%%%% SETUP CUBOID lines %%%%
lines = setup_lines_cross(particlePositions, springs, DIM);

%%%%%%%%% PARTICLE PERTRUBATIONS %%%%%%%%
%%%% INITIAL POSITION %%%%
particlePositionsOld = particlePositions;
particlePositions(:,:) = particlePositions(:,:) + [0,0,0];
%%%% INITIAL VELOCITIES %%%%
particleVelocities = zeros(particleNumbers,3);
particleVelocitiesOld = particleVelocities;
particleVelocities(:,:) = particleVelocities(:,:) + [0,0,0];
%%%% INITIAL FORCES %%%%
particleForces = zeros(particleNumbers, 3);

%%%%%%%%%%%%%%% GRAPHICS SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
%figure(1);
if DIM == 2
    axis equal; 
elseif DIM == 3
    daspect([1,1,1]);
    axis([particlePositions(1,1), particlePositions(end,1),...
        particlePositions(1,2),particlePositions(end,2),...
        -5,10]);
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
end

%%%% INITIAL ENERGY UPDATE %%%%
particleEnergyKinetic(1) = energy_kinetic(particleVelocities, particleMasses);
particleEnergySprings(1) = energy_springs(particlePositions,springs,springKS);
particleEnergyPotential(1) = energy_gravity(particlePositions,particleMasses,G,DIM);
particleEnergyTotal(1) = particleEnergyKinetic(1) + particleEnergySprings(1) + particleEnergyPotential(1);

ballEnergyKinetic(1) = energy_kinetic(ballVelocities, ballMasses);
ballEnergyPotential(1) = energy_gravity(ballPositions,ballMasses,G,DIM);
ballEnergyTotal(1) = ballEnergyKinetic(1) + ballEnergyPotential(1);

systemEnergyTotal(1) = particleEnergyKinetic(1) + ballEnergyKinetic(1);

keyboard;

%%%% INITIAL UPDATE %%%%
[particleVelocities,particleForces] = init_update(dt,particlePositions,particleVelocities,G,particleMasses,springKS,springKD,springs,DIM);
[ballVelocities,ballForces] = init_update_balls(dt,ballPositions,ballVelocities,G,ballMasses,DIM);

%profile on;

%%%% MAIN LOOP %%%%
atEnergies = 2;
updateNumber = 1;
for t = TIMES
    %half euler step forward to calculate forces in sync
    %particleVelocitiesHalfForward = particleVelocities + 0.5*dt*particleForces./particleMasses;
    %ballVelocitiesHalfForward = ballVelocities + 0.5*dt*ballForces./ballMasses;
    
    particleForces = zeros(particleNumbers,3);
    %UPDATE GRAVITY
    particleForces = calc_gravity_forces(particleForces,particleMasses,G,DIM);
    ballForces = calc_gravity_forces(ballForces, ballMasses,G,DIM);
    
    %UPDATE SPRING FORCES, use half euler forward
    particleForces = calc_spring_forces(particlePositions,particleVelocities,particleForces,springKS,springKD,springs);
    
    %UPDATE POSITION AND VELOCITY
    %particles
    [particlePositions(particleInnerIndices,:), particleVelocities(particleInnerIndices,:)] = update_rv(dt,particlePositions(particleInnerIndices,:),...
        particleVelocities(particleInnerIndices,:),particleForces(particleInnerIndices,:),particleMasses(particleInnerIndices,:));
    %balls
    [ballPositions, ballVelocities] = update_rv(dt,ballPositions,...
        ballVelocities,ballForces,ballMasses);
    
    
    %UPDATE COLLISION
    [particlePositions(particleInnerIndices,:),particleVelocities(particleInnerIndices,:),ballPositions,ballVelocities] = ...
        update_coll_mb(particlePositions(particleInnerIndices,:),particleVelocities(particleInnerIndices,:),particleForces(particleInnerIndices,:),...
        particleMasses(particleInnerIndices,:),ballPositions,ballVelocities,ballMasses,ballRadiuses);
    
    
    % take a half euler step forward to calculate energies in sync
    particleVelocitiesHalfForward = particleVelocities + 0.5*dt*particleForces./particleMasses;
    ballVelocitiesHalfForward = ballVelocities + 0.5*dt*ballForces./ballMasses;
    
    if mod(updateNumber,drawUpdate) == 0
        lines = draw_lines(particlePositions, lines);
        balls = draw_balls(balls,ballPositions,ballRadiuses);
        drawnow;
    end
    
    %%%% CALC ENERGIES %%%%
    if mod(updateNumber,energyUpdate) == 0
        energyTimes(atEnergies) = t;
        % PARTICLE ENERGIES
        particleEnergyKinetic(atEnergies) = energy_kinetic(particleVelocitiesHalfForward, particleMasses);
        particleEnergySprings(atEnergies) = energy_springs(particlePositions,springs,springKS);
        particleEnergyPotential(atEnergies) = energy_gravity(particlePositions,particleMasses,G,DIM);
        particleEnergyTotal(atEnergies) = particleEnergyKinetic(atEnergies)...
            + particleEnergySprings(atEnergies)...
            + particleEnergyPotential(atEnergies);
        
        % BALL ENERGIES
        ballEnergyKinetic(atEnergies) = energy_kinetic(ballVelocitiesHalfForward, ballMasses);
        ballEnergyPotential(atEnergies) = energy_gravity(ballPositions,ballMasses,G,DIM);
        ballEnergyTotal(atEnergies) = ballEnergyKinetic(1) + ballEnergyPotential(1);
        
        % TOTAL SYSTEM ENERGIES
        systemEnergyTotal(atEnergies) = particleEnergyTotal(atEnergies) + ballEnergyTotal(atEnergies);
        
        atEnergies = atEnergies + 1;
    end
    
    updateNumber = updateNumber + 1;
end

%profile viewer;

figure(2);
subplot(2,2,1);
plot(energyTimes, particleEnergyKinetic, '-r');
hold on;
plot(energyTimes, particleEnergySprings, '-b');
hold on;
plot(energyTimes, particleEnergyPotential, '-g');
hold on;
plot(energyTimes, particleEnergyTotal, '-m');
legend('Kinetic','springs','Potential','Total');

subplot(2,2,2);
plot(energyTimes, ballEnergyKinetic, '-r');
hold on;
plot(energyTimes, ballEnergyPotential, '-b');
hold on;
plot(energyTimes, ballEnergyTotal, '-m');
legend('Kinetic', 'Potential', 'Total');

subplot(2,1,2);
plot(energyTimes, systemEnergyTotal, '-m');
legend('System Total');

maxEnergyDifference = max_energy_diff(particleEnergyTotal)



