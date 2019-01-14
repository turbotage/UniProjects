function [R_n,V_n] = update(t,R_n,V_n, springs, M)
%UPDATE Summary of this function goes here
%   Detailed explanation goes here
    F_n = zeroes(length(R_n));
    for i=1:length(springs)
        r = R_n(springs(i).fromPI,:) - R_n(springs(i).toPI,:);
        d = norm(r);
        v = V_n(springs(i).fromPI,:) - V_n(springs(i).toPI,:);
        f = ((spring(i).ks*(d-springs(i).L) + springs(i).kd*(v.*r)/d)/d)*r; % f_ba not f_ab
        F_n(spring(i).toPI) = F_n(spring(i).toPI) + f;
        F_n(spring(i).fromPI) = F_n(spring(i).fromPI) - f;
    end
    
    V_n = V_n + F_n./M;
    
    
    
    %%TODO add gravity
    
    
    
end

