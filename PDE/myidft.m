function y = myidft(z)
% Compute the IDFT of a vector y of length N
% z k = 1/N sum {1=0}ˆ{N-1}y 1 exp(-2 pi k 1/N)
    N = length(z);
    y = zeros(1,N);
    omega = exp(2.*pi.*1i./N);
    for n=1:N
        for j=1:N
            y(n) = y(n) + z(j)*omega^((n-1)*(j-1));
        end
    end
end

