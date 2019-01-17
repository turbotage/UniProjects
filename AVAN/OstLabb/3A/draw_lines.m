function lines = draw_lines(R_n,lines)
%DRAW Summary of this function goes here
%   Detailed explanation goes here
    for i = 1:length(lines)
       set(lines(i).line, ... 
           'XData', [R_n(lines(i).fromPI, 1), R_n(lines(i).toPI, 1)],...
           'YData',[R_n(lines(i).fromPI, 2), R_n(lines(i).toPI, 2)],...
           'ZData',[R_n(lines(i).fromPI, 3), R_n(lines(i).toPI, 3)] );
    end
end

