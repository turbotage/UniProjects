
SP = [0,0,0];
NC = 2;
NR = 2;
NL = 1;
NS = NL*NC*(NR-1)+NR*NC*(NL-1)+NL*NR*(NC-1)+... %straights
     2*NL*(NR-1)*(NC-1)+2*NR*(NC-1)*(NL-1)+2*NC*(NL-1)*(NR-1)+... %2cross
     4*(NL-1)*(NR-1)*(NC-1); %3cross
Ks = 100;
Kd = 0.5;
dim = 3;

%SETUP POINTS AND SPRINGS [NC,NR,NL,S,SP,KS,KD]
[PR_n, springs] = setup_cuboid(NC,NR,NL,1,SP,Ks,Kd);

%SETUP GRAPHICS
axis([-3,3,-3,3,-3,3]);
lines = setup_lines_cross(PR_n, springs, dim);