function [p, ptot] = momentum(m, vx, vy)
    p = vx;
    ptot = zeros(length(vx(:,1)),1);
    for i = 1:length(vx(:,1))
        p(i,:) = m.*sqrt((vx(i,:).*vx(i,:) + vy(i,:).*vy(i,:)));
        ptot(i) = sum(p(i,:));
    end
    
end

