function [PR_n, PV_n] = update_coll(PR_n, PV_n, ground_balls, BR_n, Kf)
%UPDATE_COLL Summary of this function goes here
%   Detailed explanation goes here
    for i=1:length(ground_balls)
        for j=1:length(PR_n(:,1))
           r = PR_n(j,:) - BR_n(i,:);
           d = norm(r);
           %disp('Hello');
           if  d <= ground_balls(i).r
               normal = dot(PV_n(j,:),r)*r;
               % TRANSLATE
               PR_n(j,:) = PR_n(j,:) + 2*((ground_balls(i).r-d)/d)*r;
               % NEW VELOCITY
               PV_n(j,:) = PV_n(j,:) - 2*normal;
           end
        end
    end
end

