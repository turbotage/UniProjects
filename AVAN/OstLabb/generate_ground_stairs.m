function [BR_n, r_n] = generate_ground_stairs(XN, YN,N,S,SP)
    BR_n = zeros(N*XN,3);
    r_n = zeros(N*XN,1);
    for i = 0:(N-1)
        for j=1:XN
            BR_n(XN*i + j,1) = SP(1) + XN*i*S + j*S+ 0.3*rand();
            BR_n(XN*i + j,2) = SP(2) - YN*i;
            BR_n(XN*i + j,3) = SP(3);
            r_n(XN*i + j) = 0.8; %+ rand()
        end
    end
end

