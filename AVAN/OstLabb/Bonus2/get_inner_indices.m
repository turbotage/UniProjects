function innerIndices = get_inner_indices(NC,NR,NL)
    atIndex = 1;
    innerIndices = zeros(0);
    if NL == 1
        for i = 1:(NR-2)
            temp = [2:1:(NC-1)]+i*NC;
            innerIndices = [innerIndices, temp];
        end
    else
        for i = 1:(NL-2)
           for j = 1:(NR-2)
              temp = [2:1:(NC-1)]+i*NC*NR + j*NC;
              innerIndices = [innerIndices,temp];
           end
        end
    end
end

