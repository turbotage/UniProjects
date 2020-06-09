function [p, ptot] = momentum(m, vx, vy)
    p = vx;
    ptot = zeros(length(vx(:,1)),1);
    
    for i = 1:length(vx(:,1))
        p(i,:) = m.*sqrt((vx(i,:).*vx(i,:) + vy(i,:).*vy(i,:)));
    end
    
    px = sum(m.*vx,2);
    py = sum(m.*vy,2);
    
    ptot = sqrt(px.*px + py.*py);
end

