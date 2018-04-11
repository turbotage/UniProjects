function [ax, ay] = acceleration(G,m,x,y)
	ax = zeros(1,length(x));
	ay = zeros(1,length(x));
	
	acx = @(xd,yd) xd/((sqrt(xd*xd + yd*yd))^3);
	acy = @(xd,yd) yd/((sqrt(xd*xd + yd*yd))^3);
	
	for i = 1:length(m)
		for j = 1:length(m)
			if(j ~= i)
				xdiff = x(j) - x(i);
				ydiff = y(j) - y(i);
			
				ax(i) = ax(i) + G*m(j)*acx(xdiff,ydiff);
				ay(i) = ay(i) + G*m(j)*acy(xdiff,ydiff);
			end
		end
	end
end