function particles = setup_particles(particlePositions, particleRadius)
    [bx,by,bz] = sphere(10);
    bx = bx*particleRadius;
    by = by*particleRadius;
    bz = bz*particleRadius;
    for i = 1:length(particlePositions(:,1))
       particles(i) = surf(bx + particlePositions(i,1), ...
           by + particlePositions(i,2),...
           bz + particlePositions(i,3), 'FaceColor', 'r');
       hold on;
    end
end

