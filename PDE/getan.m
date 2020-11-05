function a = getan(z)
%GETAN Function that returns the values for the a_n coefficients
    N = length(z)/2;
    a = zeros(1,N);
    
    an = z(2:end)+flip(z(2:end));
    a(1) = z(1) * 2;
    a(2:end) = an(1:floor((length(an)/2)));
    a = real(a);
end