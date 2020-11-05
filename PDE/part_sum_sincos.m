function y = part_sum_sincos(z,x,N)
    [a0,a,b] = myfouriercoeff(z);
    y = ones(1,length(x))*a0/2;
    for i = 1:N
        y = y + a(i).*cos(i.*x) + b(i).*sin(i.*x);
    end
    
end

