distanceSensors = 0.1; %m

bulletMass = 0.49 * 0.001; % g

ballMass = [42.07+0.49,42.07+0.49*2,40.45+0.49,41.68+0.49,42.07+0.49*5] .* 0.001; % g
deltaT = [1.02, 0.974, 0.912, 0.968, 0.980] .* 0.001; % ms
heights = [5.9,5.7,8.6,6.8,6.1] .* 0.01; % cm

v_t = distanceSensors ./ deltaT



v_p = sqrt(2 .* ballMass .* 9.82 .* heights ./ bulletMass)

figure(1);
plot(1:1:length(deltaT), v_t);
hold on;
plot(1:1:length(deltaT), v_p);

figure(2);
plot(1:1:length(deltaT), v_p ./ v_t);

figure(3);
energyLoss = (v_t .* v_t .* bulletMass .* 0.5) - (ballMass .* 9.82 .* heights);
plot(1:1:length(deltaT),energyLoss ./ (v_t .* v_t .* bulletMass .* 0.5));
