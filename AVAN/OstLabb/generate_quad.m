function positions = generate_quad(NR,NC,S,SP)
%GENERATE_QUAD Summary of this function goes here
%   Detailed explanation goes here
    positions = zeros(NR*NC, 3);
    for i=0:(NR-1)
        for j=1:NC
            positions(i*NC+j,1) = j*S + SP(1);
            positions(i*NC+j,2) = (NR-i)*S + SP(2);
            positions(i*NC+j,3) = SP(3); %randi(5) + T(3);
        end
    end
end

