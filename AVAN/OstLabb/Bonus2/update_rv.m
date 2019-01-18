function [R_n,V_n] = update_rv(dt,R_n,V_n,F_n,M)
    V_n = V_n + dt*F_n./M;
    R_n = R_n + dt*V_n;
end

