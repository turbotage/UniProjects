function lines = draw_lines(R_n,lines)
%DRAW Summary of this function goes here
%   Detailed explanation goes here
    set(lines(1),...
        'XData',[R_n(1,1),R_n(2,1),R_n(3,1),R_n(4,1),R_n(1,1),R_n(5,1),R_n(6,1),R_n(7,1),R_n(8,1),R_n(5,1)],...
        'YData',[R_n(1,2),R_n(2,2),R_n(3,2),R_n(4,2),R_n(1,2),R_n(5,2),R_n(6,2),R_n(7,2),R_n(8,2),R_n(5,2)],...
        'ZData',[R_n(1,3),R_n(2,3),R_n(3,3),R_n(4,3),R_n(1,3),R_n(5,3),R_n(6,3),R_n(7,3),R_n(8,3),R_n(5,3)]);
    set(lines(2),...
        'XData',[R_n(2,1),R_n(6,1)],...
        'YData',[R_n(2,2),R_n(6,2)],...
        'ZData',[R_n(2,3),R_n(6,3)]);
    set(lines(3),...
        'XData',[R_n(3,1),R_n(7,1)],...
        'YData',[R_n(3,2),R_n(7,2)],...
        'ZData',[R_n(3,3),R_n(7,3)]);
    set(lines(4),...
        'XData',[R_n(4,1),R_n(8,1)],...
        'YData',[R_n(4,2),R_n(8,2)],...
        'ZData',[R_n(4,3),R_n(8,3)]);
end

