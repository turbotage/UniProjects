function z = myifft(y,N)
    y = conj(y);
    z = myfft(y,N);
    z = conj(z);
    z = z .* N;
end

