function pi = pi_map(NC,NR,NL,POS)
%PI_MAP Summary of this function goes here
%   Detailed explanation goes here
    if (POS(1) >= NC) || (POS(1) < 0) %if bigger than xout
        pi = -1;
    elseif (POS(2) >= NR) || (POS(2) < 0)
        pi = -1;
    elseif (POS(3) >= NL) || (POS(3) < 0)
        pi = -1;
    else
        % Z-part,   Y-part,   X-part
        pi = NC*NR*POS(3) + NC*POS(2) + POS(1) + 1;
    end
end

