function [x,y,vx,vy,ax,ay,t] = orbit_Nbody(G, m, x0, y0, vx0, vy0, dt, tmax)
	t = dt:dt:tmax;
    t = t - dt;
	steps = floor(tmax/dt);
    
	x = zeros(steps,length(m));
	y = zeros(steps,length(m));
	vx = zeros(steps,length(m));
	vy = zeros(steps,length(m));
	ax = zeros(steps,length(m));
	ay = zeros(steps,length(m));
	
	x(1,:) = x0;
	y(1,:) = y0;
	vx(1,:) = vx0;
	vy(1,:) = vy0;
	[atx, aty] = acceleration(G,m,x0,y0);
	ax(1,:) = atx;
	ay(1,:) = aty;
	
	for i = 1:(steps-1)
		%position
		[xt, yt] = position(x(i,:),y(i,:),vx(i,:),vy(i,:),ax(i,:),ay(i,:),dt);
		x(i+1,:) = xt;
		y(i+1,:) = yt;
		
		%acceleration
		[atx, aty] = acceleration(G,m,x(i,:),y(i,:));
		ax(i+1,:) = atx;
		ay(i+1,:) = aty;
		
		%velocity
		[vtx, vty] = velocity(vx(i,:),vy(i,:),ax(i,:),ay(i,:),ax(i+1,:),ay(i+1,:),dt);
		vx(i+1,:) = vtx;
		vy(i+1,:) = vty;
	end

end