function particles = draw_particles(particles, particlePositions, particleRadius)
    [bx,by,bz] = sphere(10);
    bx = bx*particleRadius;
    by = by*particleRadius;
    bz = bz*particleRadius;
    for i = 1:length(particlePositions(:,1))
        set(particles(i), 'Xdata', bx + particlePositions(i,1), ...
            'Ydata', by + particlePositions(i,2),...
            'Zdata', bz + particlePositions(i,3));
    end
end