function BR_n = generate_ground(NX,NY,S,SP,BR,DIM)
    if DIM == 2
        BR_n = zeros(NX,3);
        %%%% BALL POSITIONS
        for i=0:(NX-1)
            bi = i+1;
            BR_n(bi,1) = i;
            BR_n(bi,2) = 0;
            BR_n(bi,3) = 0;
        end
        BR_n = S*BR_n + SP;
        %%%% BALL DRAWING
        for i=1:length(BR_n(:,1))
            rectangle('Position', [BR_n(i,1)-BR, BR_n(i,2)-BR, ...
               BR*2, BR*2], ...
               'Curvature', [1,1], 'EdgeColor', 'r');% 'FaceColor', 'r');
        end
    elseif DIM == 3
        BR_n = zeros(NX*NY,3);
        %%%% BALL POSITIONS
        for i=0:(NY-1) %y
            for j=0:(NX-1) %x
                bi = NX*i+j+1;
                BR_n(bi,1) = j;
            	BR_n(bi,2) = i;
                BR_n(bi,3) = 0;
            end
        end
        BR_n = S*BR_n + SP;
        
        %%%% SPHERE DRAWING
        [x,y,z] = sphere(5);
        x = x*BR;
        y = y*BR;
        z = z*BR;
        for i=1:length(BR_n(:,1))
            surf(x + BR_n(i,1), y + BR_n(i,2),z + BR_n(i,3),'FaceColor','r');
            hold on;
        end
    end
    
end

