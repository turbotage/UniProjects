function [vx, vy] = velocity(vx0,vy0,ax0,ay0,ax1,ay1,dt)
	%vx = vx0;
	%vy = vy0;
    vx = vx0 + (ax0 + ax1)*dt*0.5;
    vy = vy0 + (ay0 + ay1)*dt*0.5;
	%for i = 1:length(vx0)
	%	vx(i) = vx0(i) + (ax0(i) + ax1(i))*dt*0.5;
	%	vy(i) = vy0(i) + (ay0(i) + ay1(i))*dt*0.5;
	%end
end