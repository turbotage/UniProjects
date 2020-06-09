function [x,y,vx,vy,ax,ay,t] = orbit_1body(G,M,x0,y0,vx0,vy0,dt,tmax)

    t = dt:dt:tmax;
    t = t - dt;
	steps = floor(tmax/dt);
    
    % allocate
	x = zeros(steps,1);
	y = zeros(steps,1);
	vx = zeros(steps,1);
	vy = zeros(steps,1);
	ax = zeros(steps,1);
	ay = zeros(steps,1);
    
    % conditions at start
	x(1) = x0;
	y(1) = y0;
	vx(1) = vx0;
	vy(1) = vy0;
    
    acx = @(xd,yd) xd/((sqrt(xd*xd + yd*yd))^3);
	acy = @(xd,yd) yd/((sqrt(xd*xd + yd*yd))^3);
    
	ax(1) = G*M*acx(x(1),y(1));
	ay(1) = G*M*acy(x(1),y(1));
	
	for i = 1:(steps-1)
		%position
		x(i+1) = x(i) + vx(i)*dt + ax(i)*dt*dt*0.5;
		y(i+1) = y(i) + vy(i)*dt + ay(i)*dt*dt*0.5;
		
		%acceleration
		ax(i+1) = -G*M*acx(x(i),y(i));
		ay(i+1) = -G*M*acy(x(i),y(i));
		
		%velocity
		vx(i+1) = vx(i) + (ax(i) + ax(i+1))*dt*0.5;
		vy(i+1) = vy(i) + (ay(i) + ay(i+1))*dt*0.5;
	end

end

