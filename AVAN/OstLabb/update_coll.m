function [PR_n, PV_n] = update_coll(PR_n, PV_n, BR_n, BR)
%UPDATE_COLL Summary of this function goes here
%   Detailed explanation goes here
    % OLD SLOW COLLISION CHECK
%     for i=1:length(BR_n(:,1))
%         for j=1:length(PR_n(:,1))
%            r = PR_n(j,:) - BR_n(i,:);
%            d = norm(r);
%            if  d < BR
%                normal = dot(PV_n(j,:),r)*r/(d*d);
%                % TRANSLATE
%                PR_n(j,:) = PR_n(j,:) + 2*((BR-d)/d)*r;
%                % NEW VELOCITY
%                PV_n(j,:) = PV_n(j,:) - 2*normal;
%            end
%         end
%     end
      
      % FASTER SHIT (alot faster for big systems)

%     Rd = PR_n(:,1);
%     R = PR_n;
%     for i=1:length(BR_n(:,1))
%         R = ((PR_n - BR_n(i,:)).^2).';
%         Rd = sqrt(sum(R));
%         collided = find(Rd < BR);
%         for j=1:length(collided)
%             pi = collided(j);
%             r = PR_n(pi,:) - BR_n(i,:);
%             d = norm(r);
%             normal = dot(PV_n(pi,:),r)*r/(d*d);
%             % TRANSLATE
%             PR_n(pi,:) = PR_n(pi,:) + 2*((BR-d)/d)*r;
%             % NEW VELOCITY
%             PV_n(pi,:) = PV_n(pi,:) - 2*normal;
%         end
%     end
    
    % WITH PERFECT TRANSLATION TO SURFACE OF GROUND

    Rd = PR_n(:,1);
    R = PR_n;
    for i=1:length(BR_n(:,1))
        R = ((PR_n - BR_n(i,:)).^2).';
        Rd = sqrt(sum(R));
        collided = find(Rd < BR);
        for j=1:length(collided)
            p_i = collided(j);
            r = PR_n(p_i,:) - BR_n(i,:);
            dr = norm(r);
            v = PV_n(p_i,:);
            dv = norm(v);
            
            % TRANSLATE
            c = acos(dot(-v,r)/(dr*dv));
            if c < 10^(-5)
                PR_n(p_i,:) = PR_n(p_i,:) - (BR-dr)*v/dv;
            else
                beta = asin((sin(c)*dr/BR));
                h = sin(c-beta)*BR/sin(c);
                PR_n(p_i,:) = PR_n(p_i,:) - h*v/dv;
            end
            
            % VELOCITIES
            r = PR_n(p_i,:) - BR_n(i,:);
            dr = norm(r);
            normal = (dot(-v,r)/(dr*dr))*r;
            PV_n(p_i,:) = PV_n(p_i,:) + 2*normal;
        end
    end

end

