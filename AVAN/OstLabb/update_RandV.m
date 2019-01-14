function [R_n,V_n] = update_RandV(dt,R_n,V_n, springs, M)
%UPDATE Summary of this function goes here
%   Detailed explanation goes here
    F_n = zeros(length(R_n),3);
    for i=1:length(springs)
        r = R_n(springs(i).fromPI,:) - R_n(springs(i).toPI,:);
        d = norm(r);
        v = V_n(springs(i).fromPI,:) - V_n(springs(i).toPI,:);
        f = ( (springs(i).ks*(d-springs(i).L) + springs(i).kd*dot(v,r)/d)/d)*r; % f_ba not f_ab
        F_n(springs(i).toPI,:) = F_n(springs(i).toPI,:) + f;
        F_n(springs(i).fromPI,:) = F_n(springs(i).fromPI,:) - f;
    end
    
    V_n = V_n + dt*F_n./M;
    
    R_n = R_n + dt*V_n;
    
    %%TODO add gravity
end

