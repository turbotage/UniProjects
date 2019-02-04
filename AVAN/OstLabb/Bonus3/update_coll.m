function [PR_n, PV_n,pressure] = update_coll(dt,PR_n,PV_n,PR,PM,outerCornerPos,outerCornerVel)
    pressure = 0;

    for i=1:length(PR_n(:,1))
        %x
        if (outerCornerPos(1) - PR_n(i,1)) < PR
            L = PR - (outerCornerPos(1) - PR_n(i,1));
            PR_n(i,1) = PR_n(i,1) - L;
            PV_n(i,1) = -PV_n(i,1) + outerCornerVel(1);
            pressure = pressure + abs(2*PV_n(i,1))*PM/dt;
        elseif (PR_n(i,1) - PR) < 0
            L = (PR - PR_n(i,1));
            PR_n(i,1) = PR_n(i,1) + L;
            PV_n(i,1) = -PV_n(i,1);
            pressure = pressure + abs(2*PV_n(i,1))*PM/dt;
        end
        %y
        if (outerCornerPos(2) - PR_n(i,2)) < PR
            L = PR - (outerCornerPos(2) - PR_n(i,2));
            PR_n(i,2) = PR_n(i,2) - L;
            PV_n(i,2) = -PV_n(i,2) + outerCornerVel(2);
            pressure = pressure + abs(2*PV_n(i,2))*PM/dt;
        elseif (PR_n(i,2) - PR) < 0
            L = (PR - PR_n(i,2));
            PR_n(i,2) = PR_n(i,2) + L;
            PV_n(i,2) = -PV_n(i,2);
            pressure = pressure + abs(2*PV_n(i,2))*PM/dt;
        end
        %z
        if (outerCornerPos(3) - PR_n(i,3)) < PR
            L = PR - (outerCornerPos(3) - PR_n(i,3));
            PR_n(i,3) = PR_n(i,3) - L;
            PV_n(i,3) = -PV_n(i,3) + outerCornerVel(3);
            pressure = pressure + abs(2*PV_n(i,3))*PM/dt;
        elseif (PR_n(i,3) - PR) < 0
            L = (PR - PR_n(i,3));
            PR_n(i,3) = PR_n(i,3) + L;
            PV_n(i,3) = -PV_n(i,3);
            pressure = pressure + abs(2*PV_n(i,3))*PM/dt;
        end
    end
    %particle collisions
    for i = 1:length(PR_n(:,1))
        for j = 1:length(PR_n(:,1))
            if i == j
               continue;
            end
            r_rel = PR_n(i,:) - PR_n(j,:);
            dr = norm(r_rel);
            if dr < 2*PR
                l = (2*PR - dr);
                v_rel = PV_n(i,:) - PV_n(j,:);
                PR_n(i,:) = PR_n(i,:) + l*(r_rel/dr);
                PR_n(j,:) = PR_n(j,:) - l*(r_rel/dr);
                
                PV_n(i,:) = PV_n(i,:) - (dot(v_rel,r_rel)/(dr*dr))*r_rel;
                PV_n(j,:) = PV_n(j,:) + (dot(v_rel,r_rel)/(dr*dr))*r_rel;
            end
        end
    end
    
%     %check particle particle collisions
%     for i=1:length(PR_n(:,1))
%         R = PR_n(i,:) - PR_n;
%         R(i,:) = [];
%         L = sqrt(sum(R.*R,2));
%         collided = find(L < 2*PR);
%         if ~isempty(collided)
%             v1 = PV_n(i,:);
%             v2 = PV_n(collided(1),:);
% 
%             v_rel = v1 - v2;
%             r_rel = PR_n(i,:) - PR_n(collided(1),:);
%             dr = norm(r_rel);
%             n = r_rel/dr;
%             l = (2*PR - dr);
%             
%             PR_n(i,:) = PR_n(i,:) + l*n;
%             PR_n(collided(1),:) = PR_n(collided(1),:) - l*n;
%             
%             %ignore mass terms since masses are identical
%             PV_n(i,:) = PV_n(i,:) - (dot(v_rel,r_rel)/(dr*dr))*r_rel;
%             PV_n(collided(1),:) = PV_n(collided(1),:) + (dot(v_rel,r_rel)/(dr*dr))*r_rel;
% 
%             %let them bounce out
%             %PR_n(i,:) = PR_n(i,:) + PV_n(i,:)*dt;
%             %PR_n(collided(1),:) = PR_n(collided(1),:) + PV_n(collided(1),:)*dt;
%             
%         end
%         
%     end
    
    %test wall bounces
%     L = (PR_n - [PR,PR,PR]);
%     D = L < [0,0,0];
%     collided = find(sum(D,2));
%     if ~isempty(collided)
%         PR_n(collided,:) = PR_n(collided,:) - 2*D(collided,:).*L(collided,:);
%         PV_n(collided,:) = PV_n(collided,:) - 2*D(collided,:).*PV_n(collided,:);
%     end
%     L = ((PR_n + [PR,PR,PR]) - outerCornerPos(:,:));
%     D = L > [0,0,0];
%     collided = find(sum(D,2));
%     if ~isempty(collided)
%         PR_n(collided,:) = PR_n(collided,:) - 2*D(collided,:).*L(collided,:);
%         PV_n(collided,:) = PV_n(collided,:) - 2*D(collided,:).*PV_n(collided,:);
%     end
    
end

