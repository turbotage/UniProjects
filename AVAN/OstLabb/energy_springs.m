function Es = energy_springs(PR_n,springs,KS)
%SENERGY_SPRINGS Summary of this function goes here
%   Detailed explanation goes here
    Es = 0;
    for i = 1:length(springs)
        Es = Es + (norm(PR_n(springs(i).fromPI,:) - PR_n(springs(i).toPI,:))-springs(i).L)^2;
    end
    Es = Es*0.5*KS;
end

