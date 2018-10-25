function z = myfft(y, N)
    % compute the FFT of a vector y of length N = 2m
    if N==1
        z=y;%Fourier transform of a scaler is itself
    else
        ze = myfft(y(1:2:N), N/2);% even entries %u
        zo = myfft(y(2:2:N), N/2);% odd entries  %v
        omega = ones(1,N/2);
        omega = omega.*exp(-1i.*2.*pi./N);% here i is the imaginary unit
        k = 0:((N/2)-1);
        z(1:(N/2)) = (ze + (omega.^k).*zo)./2;
        z((N/2+1):N) = (ze - (omega.^k).*zo)./2;
    end
end
