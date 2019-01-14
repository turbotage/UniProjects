function F_n = update_g(F_n,M,G,dim)
%UPDATE_G Summary of this function goes here
%   Detailed explanation goes here
    if dim == 2
        F_n(:,2) = G*M;
    elseif dim==3
        F_n(:,3) = G*M;
    end
end

