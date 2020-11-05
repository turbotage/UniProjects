function [a0,a,b] = myfouriercoeff(z)
% GET fourier coeffs  
 an = getan(z);
 a0 = an(1);
 a = an(2:end);
 b = getbn(z);
end

