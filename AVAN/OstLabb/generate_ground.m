function [BR_n, r_n] = generate_ground(N,S,SP)
%GENERATE_GROUND Summary of this function goes here
%   Detailed explanation goes here
    BR_n = zeros(N,3);
    r_n = zeros(N,1);
    for i=1:N
        BR_n(i,1) = SP(1) + (i-1)*S;
        BR_n(i,2) = SP(2);
        BR_n(i,3) = SP(3);
        r_n(i) = 0.8; %+ rand()
    end
    
end

