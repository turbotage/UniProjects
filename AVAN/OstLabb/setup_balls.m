function balls = setup_balls(BR_n,r_n, dim)
%SETUP_BALLS Summary of this function goes here
%   Detailed explanation goes here
    balls = struct([]);
    if dim == 2
        for i = 1:length(BR_n)
           balls(i).posI = i;
           balls(i).r = r_n(i);
           balls(i).ball = rectangle('Position', [BR_n(i,1)-r_n(i), BR_n(i,2)-r_n(i), ...
               r_n(i)*2, r_n(i)*2], ...
               'Curvature', [1,1], 'EdgeColor', 'r', 'FaceColor', 'w');
        end
    elseif dim == 3
        
        for i = 1:length(BR_n)
            balls(i).pos = i;
            balls(i).r = r_n(i);
            
            
            
        end
    end
end



