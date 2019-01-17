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
    
    R = BR_n;
    for i=1:length(PR_n(:,1))
        R = PR_n(i,:) - BR_n;
        collided = find(sum(R.*R,2) < BR*BR);
        if isempty(collided) == 1
           continue; 
        end
        
        r = PR_n(i,:) - BR_n(collided(1),:);
        dr = norm(r);
        v = PV_n(i,:);
        dv = norm(v);

        c = acos(dot(-v,r)/(dr*dv));
%         if isreal(c) == 0
%            disp('complex'); 
%         end
        if c < 10^(-14)
            disp('in c==0');
            % TRANSLATE
            PR_n(i,:) = PR_n(i,:) - (BR-dr)*v/dv;
            % VELOCITIES
            PV_n(i,:) = -PV_n(i,:);
            
            dt = (BR-dr)/dv;
            PR_n(i,:) = PR_n(i,:) + PV_n(i,:)*dt;
        else
            beta = asin((sin(c)*dr/BR));
            h = sin(c-beta)*BR/sin(c);
            % TRANSLATE
            PR_n(i,:) = PR_n(i,:) - h*v/dv;

            %TEST IF WE STILL ARE INSIDE A BALL
            %R = PR_n(i,:) - BR_n;
            %collided = find(sum(R.*R,2) < BR*BR)
            %seems to not be needed

            % VELOCITIES
            r = PR_n(i,:) - BR_n(collided(1),:);
            dr = norm(r);
            normal = (dot(v,r)/(dr*dr))*r;
            PV_n(i,:) = PV_n(i,:) - 2*normal;
            
            %LET THE PARTICLE BOUNCE AWAY FROM THE SURFACE
            dv = norm(PV_n(i,:));
            dt = h/dv;
            PR_n(i,:) = PR_n(i,:) + dt*PV_n(i,:);
        end
    end

end

