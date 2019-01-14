function lines = setup_lines(R_n, springs, NR, NC)
    lines = struct([]);
    % TOP
    k = 1;
    for i = 1:(NC-1)
        lines(k).fromPI = springs(i).fromPI;
        lines(k).toPI = springs(i).toPI;
        lines(k).line = line(...
            [R_n(springs(i).fromPI,1), R_n(springs(i).toPI,1)], ...
            [R_n(springs(i).fromPI,2), R_n(springs(i).toPI,2)], ...
            [R_n(springs(i).fromPI,3), R_n(springs(i).toPI,3)]);
        %set(lines(i).line,'erasemode','xor');
        k = k + 1;
    end
    % BOTTEM
    for j = 1:(NC-1)
        i = (NC-1)*(NR-1) + j;
        lines(k).fromPI = springs(i).fromPI;
        lines(k).toPI = springs(i).toPI;
        lines(k).line = line(...
            [R_n(springs(i).fromPI,1), R_n(springs(i).toPI,1)], ...
            [R_n(springs(i).fromPI,2), R_n(springs(i).toPI,2)], ...
            [R_n(springs(i).fromPI,3), R_n(springs(i).toPI,3)]);
        %set(lines(i).line,'erasemode','xor');
        k = k + 1;
    end
    % LEFT LINE
    for j = 1:(NR-1)
        i = (NC-1)*NR+(j-1)*NC + 1;
        lines(k).fromPI = springs(i).fromPI;
        lines(k).toPI = springs(i).toPI;
        lines(k).line = line(...
            [R_n(springs(i).fromPI,1), R_n(springs(i).toPI,1)], ...
            [R_n(springs(i).fromPI,2), R_n(springs(i).toPI,2)], ...
            [R_n(springs(i).fromPI,3), R_n(springs(i).toPI,3)]);
        %set(lines(i).line,'erasemode','xor');
        k = k + 1;
    end
    % RIGHT LINE
    for j = 1:(NR-1)
        i = (NC-1)*NR+j*NC;
        lines(k).fromPI = springs(i).fromPI;
        lines(k).toPI = springs(i).toPI;
        lines(k).line = line(...
            [R_n(springs(i).fromPI,1), R_n(springs(i).toPI,1)], ...
            [R_n(springs(i).fromPI,2), R_n(springs(i).toPI,2)], ...
            [R_n(springs(i).fromPI,3), R_n(springs(i).toPI,3)]);
        %set(lines(i).line,'erasemode','xor');
        k = k + 1;
    end
end