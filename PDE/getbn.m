function b = getbn(z)
%GETAN Function that returns the values for the b_n coefficients
    N = length(z)/2;
    b = zeros(1,N);
    
    bn = 1i.*(z(2:end)-flip(z(2:end)));
    %b(2:end) = bn(1:floor((length(bn)/2)));
    b = bn(1:floor(length(bn)/2));
    b = real(b);
end