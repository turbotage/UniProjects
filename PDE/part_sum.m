function y = part_sum(z, x)
    omega = exp(1i);
    y = zeros(1,length(x));
    N = length(z);
    for n = 1:N/2
        y = y + z(n).*omega.^((n-1).*x); 
    end
    for n = N/2+1:N
        y = y + z(n).*omega.^((-N+(n-1)).*x);
    end
    
end

