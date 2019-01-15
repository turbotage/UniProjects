function [PR_n,springs] = setup_cuboid(NC,NR,NL,S)
    %GENERATE POSITIONS
    PR_n = zeros(NR*NC*NL, 3);
    for k=0:(NL-1) %z
        for j=0:(NR-1) %y
            for i=0:(NC-1) %z
                pi = NC*NR*k + NC*j + i + 1;
                PR_n(pi,1) = i;
                PR_n(pi,2) = j;
                PR_n(pi,3) = k;
            end
        end
    end
    %SETUP SPRINGS INDICES
    springs = struct('fromPI',{},'toPI',{},'L',{});
    SN = 1;
    for k = 0:(NL-1) %z
        for j = 0:(NR-1) %y
            for i = 0:(NC-1) %x
                %%%% STRAIGHTS
                %x
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i+1,j,k]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                %y
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i,j+1,k]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                %z
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i,j,k+1]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                %%%% 2CROSS
                % 2cross x+ y+ z0
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i+1,j+1,k]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                % 2cross x- y+ z0
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i-1,j+1,k]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                % 2cross x+ y0 z+
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i+1,j,k+1]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                % 2cross x- y0 z+
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i-1,j,k+1]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                % 2cross x0 y+ z+
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i,j+1,k+1]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                % 2cross x0 y- z+
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i,j-1,k+1]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                %%%% 3CROSS
                % 3cross x+ y+ z+
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i+1,j+1,k+1]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                % 3cross x- y+ z+
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i-1,j+1,k+1]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                % 3cross x+ y- z+
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i+1,j-1,k+1]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
                % 3cross x- y- z+
                from = pi_map(NC,NR,NL,[i,j,k]);
                to = pi_map(NC,NR,NL,[i-1,j-1,k+1]);
                if to ~= -1
                   springs(SN) = setup_spring(from,to);
                   SN = SN+1;
                end
            end
        end
    end
    
    %MAKE TRANSLATIONS AND SCALE FOR PARTICLE POSITIONS
    PR_n = S*PR_n; % scale, translation
    
    %SETUP SPRING COEEF AND DISTANCES
    for i = 1:length(springs)
       springs(i).L = norm(PR_n(springs(i).fromPI,:) - PR_n(springs(i).toPI,:));
    end
    
end

