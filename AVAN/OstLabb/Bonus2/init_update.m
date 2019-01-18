function [V_n,F_n] = init_update(dt,R_n,V_n,G,M,KS,KD,springs,DIM)
    F_n = zeros(length(R_n(:,1)),3);
    
    if DIM == 2
        F_n(:,2) = F_n(:,2) + G.*M;
    elseif DIM==3
        F_n(:,3) = F_n(:,3) + G.*M;
    end
    
    for i=1:length(springs)
        r = R_n(springs(i).fromPI,:) - R_n(springs(i).toPI,:);
        d = norm(r);
        v = V_n(springs(i).fromPI,:) - V_n(springs(i).toPI,:);
        f = ( (KS*(d-springs(i).L) + KD*dot(v,r)/d)/d)*r; % f_ba not f_ab
        F_n(springs(i).toPI,:) = F_n(springs(i).toPI,:) + f;
        F_n(springs(i).fromPI,:) = F_n(springs(i).fromPI,:) - f;
    end
    
    V_n = V_n - 0.5*dt*F_n./M;
end

