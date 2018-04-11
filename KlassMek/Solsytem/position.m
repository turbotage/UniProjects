function [x,y] = position(x0, y0, vx, vy, ax, ay, dt)
	x = x0;
	y = y0;
	for i = 1:length(x0)
		x(i) = x0(i) + vx(i)*dt + ax(i)*dt*dt/2;
		y(i) = y0(i) + vy(i)*dt + ay(i)*dt*dt/2;
	end
end