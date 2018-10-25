function b = getbn(z)
%GETAN Summary of this function goes here
%   Detailed explanation goes here
    N = length(z)/2;
    a = zeros(1,N-1);
    
    an = 1i.*z(2:end)+flip(z(1:end-1));
    a(2:end) = an(1:floor((length(an)/2)));
end