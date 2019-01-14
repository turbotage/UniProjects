function positions = generate_quad(NR,NC,S,T)
%GENERATE_QUAD Summary of this function goes here
%   Detailed explanation goes here
    positions = ones(NR*NC, 3);
    for i=0:(NR-1)
        for j=1:NC
            positions(i*NC+j,1) = j*S + T(1);
            positions(i*NC+j,2) = (NR-i)*S + T(2);
            positions(i*NC+j,3) = randi(5) + T(3);
        end
    end
end

