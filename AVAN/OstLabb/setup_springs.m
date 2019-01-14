function springs = setup_springs(R_n,NR,NC, Ks, Kd)
%SETUP_SPRINGS Summary of this function goes here
%   Detailed explanation goes here
    springs = struct([]);
    i = 0;
    %x-axis springs
    for j=0:(NR-1)
        for k=1:(NC-1)
            i=i+1; % used so that the springs gets an own slot
            pi = j*NC + k;
            springs(i).fromPI = pi;
            springs(i).toPI = pi+1;
            springs(i).ks = Ks;
            springs(i).kd = Kd;
            springs(i).L = norm(R_n(pi) - R_n(pi+1));
            %springs(i).line = line([R_n(pi,1), R_n(pi+1,1)], ...
            %    [R_n(pi,2), R_n(pi+1,2)], [R_n(pi,3), R_n(pi+1,3)]);
        end
    end
    %y-axis springs
    for pi=1:NC*(NR-1)
        i=i+1;
        springs(i).fromPI = pi;
        springs(i).toPI = pi+NC;
        springs(i).ks = Ks;
        springs(i).kd = Kd;
        springs(i).L = norm(R_n(pi) - R_n(pi+NC));
        %springs(i).line = line([R_n(pi,1), R_n(pi+NC,1)], ...
        %    [R_n(pi,2), R_n(pi+NC,2)], [R_n(pi,3), R_n(pi+NC,3)]);
    end
    % right diagonal
    for j=0:(NR-2)
        for k=1:(NC-1)
            i=i+1;
            pi = j*NC + k;
            springs(i).fromPI = pi;
            springs(i).toPI = pi+NC+1;
            springs(i).ks = Ks;
            springs(i).kd = Kd;
            springs(i).L = norm(R_n(pi) - R_n(pi+NC+1));
            %springs(i).line = line([R_n(pi,1), R_n(pi+NC+1,1)], ...
            %    [R_n(pi,2), R_n(pi+NC+1,2)], [R_n(pi,3), R_n(pi+NC+1,3)]);
        end
    end
    % left diagonal
    for j=0:(NR-2)
        for k=2:NC
            i=i+1;
            pi = j*NC + k;
            springs(i).fromPI = pi;
            springs(i).toPI = pi+NC-1;
            springs(i).ks = Ks;
            springs(i).kd = Kd;
            springs(i).L = norm(R_n(pi) - R_n(pi+NC-1));
            %springs(i).line = line([R_n(pi,1), R_n(pi+NC-1,1)], ...
            %    [R_n(pi,2), R_n(pi+NC-1,2)], [R_n(pi,3), R_n(pi+NC-1,3)]);
        end
    end
end

