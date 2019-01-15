function [PR_n, PV_n] = update_coll(PR_n, PV_n, BR_n, BR)
%UPDATE_COLL Summary of this function goes here
%   Detailed explanation goes here
    for i=1:length(BR_n(:,1))
        for j=1:length(PR_n(:,1))
           r = PR_n(j,:) - BR_n(i,:);
           d = norm(r);
           if  d < BR
               normal = dot(PV_n(j,:),r)*r/(d*d);
               % TRANSLATE
               PR_n(j,:) = PR_n(j,:) + 2*((BR-d)/d)*r;
               % NEW VELOCITY
               PV_n(j,:) = PV_n(j,:) - 2*normal;
           end
        end
    end
end

