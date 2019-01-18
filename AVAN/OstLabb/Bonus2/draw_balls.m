function balls = draw_balls(balls, ballPositions, ballRadiuses)
    [bx,by,bz] = sphere;
    for i = 1:length(balls)
       bx = bx*ballRadiuses(i);
       by = by*ballRadiuses(i);
       bz = bz*ballRadiuses(i);
       set(balls(i), 'XData', bx + ballPositions(i,1), ...
        'YData', by + ballPositions(i,2),...
        'ZData', bz + ballPositions(i,3)); 
    end
end

