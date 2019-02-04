function [V_n,F_n] = init_update(dt,R_n,V_n,G,M,DIM)
    F_n = zeros(length(R_n(:,1)),3);
    
    if DIM == 2
        F_n(:,2) = F_n(:,2) + G*M;
    elseif DIM==3
        F_n(:,3) = F_n(:,3) + G*M;
    end
    
    V_n = V_n - 0.5*dt*F_n./M;
end

