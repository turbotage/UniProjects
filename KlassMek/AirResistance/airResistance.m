%order of csv is, radius,frequency

dimplesDiameter = 4.17 * 0.01; % cm +- 0.05 mm
dimplesRadius = dimplesDiameter * 0.5;
smoothDiameter = 3.76 * 0.01; % cm +- 0.05 mm
smoothRadius = smoothDiameter * 0.5;
axisDiameter = 2.06 * 0.01; % cm +- 0.1 mm
axisRadius = axisDiameter * 0.5;

dimplesArea = pi*dimplesRadius*dimplesRadius;
smoothArea = pi*smoothRadius*smoothRadius;

rho = 1.21; % kg/m^3
eta = 1.86*10^(-5); % Ns/m^2

smoothMass = 2.5 * 0.001; % g +- 20 mg
dimplesMass = 4.4 * 0.001; % g +- 20 mg

dimples = load('dimples.csv');
smooth = load('smooth.csv');

% subtract ball radius
dimplesR = (dimples(:,1) .* 0.01) - dimplesRadius;
smoothR = (smooth(:,1) .* 0.01) - smoothRadius;

% divide by three because values 
% in csv is motor frequency 
dimplesOmega = (dimples(:,2) ./ 3) .* (2*pi);
smoothOmega = (smooth(:,2) ./ 3) .* (2*pi);

%velocities
dimplesVelocity = dimplesR .* dimplesOmega;
smoothVelocity = smoothR .* smoothOmega;

cds_num = 2 .* smoothMass .* axisRadius;
cds_den = (rho * smoothArea) .* smoothR .* sqrt(smoothR.^2 - axisRadius.^2);
cdd_num = 2 .* dimplesMass .* axisRadius;
cdd_den = (rho * dimplesArea) .* dimplesR .* sqrt(dimplesR.^2 - axisRadius.^2);

Cd_dimples = cdd_num ./ cdd_den;
Cd_smooth = cds_num ./ cds_den;

Re_dimples = rho .* dimplesVelocity .* dimplesDiameter ./ eta;
Re_smooth = rho .* smoothVelocity .* smoothDiameter ./ eta;

figure(1);
loglog(Re_dimples, Cd_dimples);
hold on;
loglog(Re_smooth, Cd_smooth);
title('Cd as function of Re');
legend('dimples', 'smooth');
xlabel('Re');
ylabel('Cd');

figure(2);
loglog(dimplesVelocity, Cd_dimples);
hold on;
loglog(smoothVelocity, Cd_smooth);
title('Cd as function of v');
legend('dimples', 'smooth');
xlabel('v');
ylabel('Cd');

f = @(x) 518.57;

figure(3);
cdv2 = Cd_dimples.*dimplesVelocity.^2
%k = sqrt((2*10^(-3)*45.73*9.82) / (1.21 * (pi/4) * 10^(-6) * 42.69^2))
%error = dimplesVelocity - (1./sqrt(Cd_dimples)).*k;
steps = 1:1:length(dimplesVelocity);
plot(steps, dimplesR);
polyfit(steps,cdv2,2)








