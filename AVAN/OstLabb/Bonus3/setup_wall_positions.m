function wallPositions = setup_wall_positions(outerCorner)
    x = outerCorner(1);
    y = outerCorner(2);
    z = outerCorner(3);
    wallPositions = ones(8,3);
    wallPositions(1,:) = [0,0,0];
    wallPositions(2,:) = [x,0,0];
    wallPositions(3,:) = [x,y,0];
    wallPositions(4,:) = [0,y,0];
    wallPositions(5,:) = [0,0,z];
    wallPositions(6,:) = [x,0,z];
    wallPositions(7,:) = [x,y,z];
    wallPositions(8,:) = [0,y,z];
end

