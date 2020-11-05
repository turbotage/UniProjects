function error = part_sum_max_error(z,x,f,N)
    part_sum = part_sum_sincos(z,x,N);
    error = max(part_sum - f(x));
end

