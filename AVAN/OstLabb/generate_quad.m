function PR_n = generate_quad(NR,NC,S,SP)
%GENERATE_QUAD Summary of this function goes here
%   Detailed explanation goes here
    PR_n = zeros(NR*NC, 3);
    for i=0:(NR-1)
        for j=1:NC
            PR_n(i*NC+j,1) = (j-1);
            PR_n(i*NC+j,2) = i;
            PR_n(i*NC+j,3) = 0; %randi(5) + T(3);
        end
    end
    PR_n = S*PR_n;
    PR_n = PR_n + SP;
end

