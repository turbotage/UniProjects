theta = -30*pi/180;

rotation = [cos(theta), -sin(theta);
            sin(theta), cos(theta)
            ];
        
origo = [1,1];
        
points = [
    1,2;
    2,4;
    

         ]

for i = 1:length(points)
    points(i) = points(i) - origo
end
      
rotation*p1
