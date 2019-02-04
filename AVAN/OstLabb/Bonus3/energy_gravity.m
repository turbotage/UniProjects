function Eg = energy_gravity(PR_n,M,G,dim)
%ENERGY_GRAVITY Summary of this function goes here
%   Detailed explanation goes here
    if dim == 2
        Eg = -1*G*sum(PR_n(:,2).*M);
    elseif dim == 3
        Eg = -1*G*sum(PR_n(:,3).*M);
    end
end

