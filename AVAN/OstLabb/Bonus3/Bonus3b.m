clf;
clear all;
%%%%%%%%%%%% INFORMATION %%%%%%%%%%%%%%%%


%%%%%%%%%%%%% GLOBAL DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%
%%%% GENERAL SETTINGS %%%%
dt = 2*10^-13;
drawUpdate = 1; %update every drawUpdate:th frame
energyUpdate = 20; %update every energyUpdate:th frame

TIMES = 0:dt:4*10^(-10);

DIM = 3; %decides how to draw, setup ground and in which direction gravity acts
G = 0; %gravity constant

T = 293; %K
P = 101325; %Pa
k = 1.38065*10^(-23);
totalPressure = 0;

%%%% PARTICLES SETTINGS %%%%
particleNumber = 50;
particleMass = 4.8*10^(-26);
particleMasses = ones(particleNumber,1)*particleMass;
particleRadius = 4*10^(-10);


%%%% ENERGY SETTINGS %%%%
energyCalcN = floor(length(TIMES)/energyUpdate)+1;
energyTimes = zeros(energyCalcN, 1);

particleEnergyKinetic = zeros(energyCalcN, 1);
particleEnergyPotential = zeros(energyCalcN, 1);
particleEnergyTotal = zeros(energyCalcN, 1);


%%%%%%%%% WALL SETTINGS %%%%%%%%
start_scale = 4;
Vfinal = [1.2591*10^(-8), 1.2591*10^(-8), 1.2591*10^(-8)];
cornerEdge = start_scale*Vfinal;
edgeVelocity = (Vfinal-cornerEdge)/TIMES(end);

%%%%%%%%% PARTICLE PERTRUBATIONS %%%%%%%%
%%%% INITIAL POSITION %%%%
particlePositions = rand(particleNumber,3).*cornerEdge*0.75 + particleRadius*[1,1,1];
%%%% INITIAL VELOCITIES %%%%
%particleVelocities = 2*rand(particleNumber,3)-1;
%velocityNorms = vecnorm(particleVelocities.').';
%particleVelocities = particleVelocities./velocityNorms;
sigma = sqrt(0.5*2*k*T/particleMass);
%v_mean = sqrt(3*k*T/particleMass);
%c = (particleMass/(2*pi*k*T))^(3/2);
%particleVelocitiesMag = abs(normrnd(0,sigma,particleNumber,1));
%particleVelocities = particleVelocities.*particleVelocitiesMag;
particleVelocities = normrnd(0,sigma,particleNumber,3);

T = energy_kinetic(particleVelocities,particleMasses)*2/(3*k)/particleNumber;
%%%% INITIAL FORCES %%%%
particleForces = zeros(particleNumber, 3);

%%%%%%%%% GUI OBJECTS %%%%%%%
particles = setup_particles(particlePositions,particleRadius);

wallPositions = setup_wall_positions(cornerEdge);
lines = setup_lines(wallPositions);

%%%%%%%%%%%%%%% GRAPHICS SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
if DIM == 2
    axis equal; 
elseif DIM == 3
    daspect([1,1,1]);
    axis([0,cornerEdge(1),0,cornerEdge(2),0,cornerEdge(3)]);
    xlabel('x-axis');
    ylabel('y-axis');
    zlabel('z-axis');
end

%%%% INITIAL ENERGY UPDATE %%%%
particleEnergyKinetic(1) = energy_kinetic(particleVelocities, particleMasses);
particleEnergyPotential(1) = energy_gravity(particlePositions,particleMasses,G,DIM);
particleEnergyTotal(1) = particleEnergyKinetic(1) + particleEnergyPotential(1);

keyboard;

%%%% INITIAL UPDATE %%%%
[particleVelocities,particleForces] = init_update(dt,particlePositions,particleVelocities,G,particleMasses,DIM);

%profile on;

%%%% MAIN LOOP %%%%
atEnergies = 2;
updateNumber = 1;
for t = TIMES
    cornerEdge = cornerEdge + edgeVelocity*dt;
    wallPositions = setup_wall_positions(cornerEdge);
    
    %UPDATE POSITION AND VELOCITY
    %particles
    [particlePositions, particleVelocities] = update_rv(dt,particlePositions,...
        particleVelocities,particleForces,particleMasses);
    
    %UPDATE COLLISION
    [particlePositions,particleVelocities,pressure] = update_coll(dt,particlePositions,particleVelocities,...
        particleRadius,particleMass,cornerEdge,edgeVelocity);
    
    totalPressure = totalPressure + pressure;
    
    % take a half euler step forward to calculate energies in sync
    particleVelocitiesHalfForward = particleVelocities + 0.5*dt*particleForces./particleMasses;
    
    if mod(updateNumber,drawUpdate) == 0
        %lines = draw_lines(particlePositions, lines);
        draw_particles(particles,particlePositions,particleRadius);
        lines = draw_lines(wallPositions,lines);
        drawnow;
    end
    
    %%%% CALC ENERGIES %%%%
    if mod(updateNumber,energyUpdate) == 0
        energyTimes(atEnergies) = t;
        % PARTICLE ENERGIES
        particleEnergyKinetic(atEnergies) = energy_kinetic(particleVelocitiesHalfForward, particleMasses);
        particleEnergyPotential(atEnergies) = energy_gravity(particlePositions,particleMasses,G,DIM);
        particleEnergyTotal(atEnergies) = particleEnergyKinetic(atEnergies)...
            + particleEnergyPotential(atEnergies);
        
        atEnergies = atEnergies + 1;
    end
    
    updateNumber = updateNumber + 1;
end


%profile viewer;

% figure(2);
% plot(energyTimes, particleEnergyKinetic, '-r');
% hold on;
% plot(energyTimes, particleEnergyPotential, '-g');
% hold on;
% plot(energyTimes, particleEnergyTotal, '-m');
% legend('Kinetic','Potential','Total');

figure(3);
plot(energyTimes, (2/(k*3))*particleEnergyKinetic/particleNumber);
legend('Temperature (K)');


