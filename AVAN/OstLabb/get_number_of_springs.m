function number_of_springs = get_number_of_springs(NC,NR,NL)
    number_of_springs = NL*NC*(NR-1)+...
    NR*NC*(NL-1)+NL*NR*(NC-1)+... %straights
    2*NL*(NR-1)*(NC-1)+...
    2*NR*(NC-1)*(NL-1)+2*NC*(NL-1)*(NR-1)+... %2cross
    4*(NL-1)*(NR-1)*(NC-1); %3cross
end

