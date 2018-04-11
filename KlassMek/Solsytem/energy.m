function [ekt, ept, ett] = energy(G,m,x,y,vx,vy)
    ek = zeros(length(x(:,1)), length(x(1,:)));
    ep = zeros(length(x(:,1)), length(x(1,:)));
    et = zeros(length(x(:,1)), length(x(1,:)));
    
    ekt = zeros(length(x(:,1)),1);
    ept = zeros(length(x(:,1)),1);
    ett = zeros(length(x(:,1)),1);
    
    nr = length(x(:,1));
    nc = length(x(1,:));
    
    for i=1:nr
        ek(i,:) = m.*(vx(i,:).*vx(i,:) + vy(i,:).*vy(i,:)).*0.5;
    end
    
    parfor i=1:nr
        for j=1:nc
            for k=(j+1):nc
                r = sqrt((x(i,j)-x(i,k))^2 +(y(i,j)-y(i,k))^2);
                ept(i) = ept(i) + (-1*G*(m(j)*m(k))/r);
            end
        end
    end
    
    for i = 1:length(x(:,1))
        ekt(i) = sum(ek(i,:));
    end
    
    ett = ekt + ept;
end