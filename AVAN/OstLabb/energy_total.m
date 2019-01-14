function E = energy_total(R_n,V_n, M)
%ENERGY_TOTAL Summary of this function goes here
%   Detailed explanation goes 
    % KINETIC
    Ek = 0.5*sum( (V_n(:,1).^2 + V_n(:,2).^2 + V_n(:,3).^2).*M );
    
    % SPRING POTENTIAL
    Es = 0;
    for i = 1:length(springs)
        Es = Es + norm(R_n(springs(i).fromPI,:) - R_n(springs(i).toPI,:))^2;
    end
    Es = Es*0.5*springs(1).KS;
    E = Ek + Es;
    
end

