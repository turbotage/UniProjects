function [prn,pvn,brn,bvn] = update_coll_mb(prn,pvn,pfn,pm,brn,bvn,bm,bR)
%     R = prn;
%     for i=1:length(brn(:,1))
%         R = brn(i,:) - prn;
%         collided = find(sum(R.*R,2) < bR(i)*bR(i));
%         if isempty(collided) == 1
%             continue;
%         end
%         
%         times = collided;
%         Ptot = [0,0,0];
%         
%         
%         % go through particles find lengths and times
%         for j=1:length(collided)
%             p_i = collided(j);
%             r = prn(p_i,:) - brn(i,:);
%             dr = norm(r);
%             v_on_ball = pvn(p_i,:) - bvn(i,:);
%             dv_on_ball = norm(v_on_ball);
%             
%             c = acos(dot(-v_on_ball,r)/(dr*dv_on_ball));
%             
%             if c < 10^(-14)
%                disp('straight collision'); 
%             else
%                 beta = asin((sin(c)*dr/bR(i)));
%                 h = sin(c-beta)*bR(i)/sin(c);
%                 dt = h/dv_on_ball;
%                 
%                 prn(p_i,:) = prn(p_i,:) - h*(v_on_ball/dv_on_ball);
%                 
%                 r = prn(p_i,:) - brn(i,:);
%                 dr = norm(r);
%                 
%                 %first add zeroing of momentum for particle
%                 Ptot = Ptot + pvn(p_i,:)*pm(p_i);
%                 
%                 %pvn(p_i,:) = norm(pvn(p_i,:))*r/dr;
%                 pvn(p_i,:) = norm(bvn(i))*r/dr;
%                 %prn(p_i,:) = prn(p_i,:) + dt*pvn(p_i,:);
%                 
%                 %add change of momentum for particle
%                 Ptot = Ptot - pvn(p_i,:)*pm(i);
%                 
%                 %add impuls from spring forces
%                 %Ptot = Ptot + dt*dot(pfn(p_i,:),r/dr)*pfn(p_i,:)/norm(pfn(p_i,:));
%                 
%             end
%         end
%         
%         bv = Ptot/bm(i);
%         bvn(i,:) = bvn(i,:) + bv;
%         
%    
%     end

    R = prn;
    for i=1:length(brn(:,1))
        R = brn(i,:) - prn;
        collided = find(sum(R.*R,2) < bR(i)*bR(i));
        if isempty(collided) == 1
            continue;
        end
        
        times = collided;
        bv_contrib = zeros(1,3);
        Ptot = zeros(1,3);
        
        % go through particles find lengths and times
        for j=1:length(collided)
            p_i = collided(j);
            r = prn(p_i,:) - brn(i,:);
            dr = norm(r);
            v_on_ball = pvn(p_i,:) - bvn(i,:);
            dv_on_ball = norm(v_on_ball);
            
            c = acos(dot(-v_on_ball,r)/(dr*dv_on_ball));
            
            if c < 10^(-14)
               disp('straight collision'); 
            else
                beta = asin((sin(c)*dr/bR(i)));
                h = sin(c-beta)*bR(i)/sin(c);
                dt = h/dv_on_ball;
                
                prn(p_i,:) = prn(p_i,:) - h*(v_on_ball/dv_on_ball);
                
                r = prn(p_i,:) - brn(i,:);
                dr = norm(r);
                
                
                pvn(p_i,:) = pvn(p_i,:) - 2*bm*(dot(v_on_ball,r)/(dr*dr))*r/(bm+pm(p_i));
                bv_contrib = bv_contrib + 2*pm(p_i)*(dot(v_on_ball,r)/(dr*dr))*r/(bm+pm(p_i));
                
                %prn(p_i,:) = prn(p_i,:) + 0.5*dt*pvn(p_i,:);
                
                %add impuls from spring forces
                Ptot = Ptot + dt*dot(pfn(p_i,:),r/dr)*pfn(p_i,:)/norm(pfn(p_i,:));
                
            end
        end
        
        bvn(i,:) = bvn(i,:) + bv_contrib - Ptot/bm;
        
    end


end

