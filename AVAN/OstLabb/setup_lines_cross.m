function lines = setup_lines_cross(R_n, springs, NR, NC)
    lines = struct([]);
    for i = 1:length(springs)
       lines(i).fromPI = springs(i).fromPI;
       lines(i).toPI = springs(i).toPI;
       lines(i).line = line(...
            [R_n(springs(i).fromPI,1), R_n(springs(i).toPI,1)], ...
            [R_n(springs(i).fromPI,2), R_n(springs(i).toPI,2)], ...
            [R_n(springs(i).fromPI,3), R_n(springs(i).toPI,3)]);
        %set(lines(i).line,'erasemode','xor');
    end
end