function F_n = calc_gravity_forces(F_n,M,G,DIM)
%UPDATE_G Summary of this function goes here
%   Detailed explanation goes here
    if DIM == 2
        F_n(:,2) = F_n(:,2) + G*M;
    elseif DIM==3
        F_n(:,3) = F_n(:,3) + G*M;
    end
end

