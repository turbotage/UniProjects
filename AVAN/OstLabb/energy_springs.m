function Es = energy_springs(springs)
%SENERGY_SPRINGS Summary of this function goes here
%   Detailed explanation goes here
    Es = 0;
    for i = 1:length(springs)
        Es = Es + norm(R_n(springs(i).fromPI,:) - R_n(springs(i).toPI,:))^2;
    end
    Es = Es*0.5*springs(1).KS;
end

