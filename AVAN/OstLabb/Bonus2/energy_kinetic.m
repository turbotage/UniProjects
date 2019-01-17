function Ek = energy_kinetic(V_n, M)
%ENERGY_KINETIC Summary of this function goes here
%   Detailed explanation goes here
    Ek = 0.5*sum( (V_n(:,1).^2 + V_n(:,2).^2 + V_n(:,3).^2).*M );
    %Ek = 0.5*sum(sum(((V_n.^2)).*M));
end

