function b = getbn(z)
%GETAN Summary of this function goes here
%   Detailed explanation goes here
    N = length(z)/2;
    b = zeros(1,N);
    
    bn = 1i.*z(2:end)+flip(z(1:end-1));
    b(2:end) = bn(1:floor((length(bn)/2)));
end