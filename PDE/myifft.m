function z = myifft(y,N)
    % compute the FFT of a vector y of length N = 2m
    y = conj(y);
    z = myfft(y,N);
    z = conj(z);
    z = z .* N;
end

