function diff = max_energy_diff(E)
%MAX_ENERGY_DIFF Summary of this function goes here
%   Detailed explanation goes here
    max = 0;
    for i = 2:length(E)
       diff = abs(1-E(1)/E(i));
       if diff > max
           max = diff;
       end
    end
    diff = max*100;
end

