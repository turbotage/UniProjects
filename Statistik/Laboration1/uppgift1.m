N = 50;

% b)

probYgt25 = 1 - binocdf(25,N,0.4);

% c)

probYget25 = 1 - binocdf(24,N,0.4);

% d)

probYe15 = binopdf(15,N,0.4);

% c)

probYmod3a5 = binopdf(15,N,0.4) + binopdf(30,N,0.4) + binopdf(45,N,0.4);

% f)

probYmod3o5 = 0;

for i = 1:N
   if (mod(i,3) == 0)
       probYmod3o5 = probYmod3o5 + binopdf(i,N,0.4);
   elseif (mod(i,5) == 0)
       probYmod3o5 = probYmod3o5 + binopdf(i,N,0.4);
   end
end

disp(probYmod3o5)

% g)

probYg20le30 = binocdf(30,N,0.4) - binocdf(19,N,0.4)