

syms x y;
f = @(x,y) y - sin(x) - 5;
g= @(x,y) y.^2 + 12*x - 4;

fg = cell(2,1);
fg{1} = f;
fg{2} = g;

format SHORT;

xy = [x y];
val1 = -1.000;
val2 = 4.5000;

iter0 = [val1, val2]
[val1, val2] = nextIteration(fg, xy, [val1, val2]);
iter1 = [val1, val2]
[val1, val2] = nextIteration(fg, xy, [val1, val2]);
iter2 = [val1, val2]
[val1, val2] = nextIteration(fg, xy, [val1, val2]);
iter3 = [val1, val2]

function [nextX, nextY] = nextIteration(funcs,vars,poss)
	
	f = funcs{1};
	g = funcs{2};
	f1 = diff(funcs(1),vars(1));
	f2 = diff(funcs(1),vars(2));
	g1 = diff(funcs(2),vars(1));
	g2 = diff(funcs(2),vars(2));

	fg2 = (atPoint(f,vars,poss)*atPoint(g2,vars,poss));
	f2g = (atPoint(f2,vars,poss)*atPoint(g,vars,poss));
	f1g2 = (atPoint(f1,vars,poss)*atPoint(g2,vars,poss));
	f2g1 = (atPoint(f2,vars,poss)*atPoint(g1,vars,poss));
	
	nextX = poss(1) - ((fg2 - f2g)/(f1g2 - f2g1));
	
	f1g = (atPoint(f1,vars,poss)*atPoint(g,vars,poss));
	fg1 = (atPoint(f,vars,poss)*atPoint(g1,vars,poss));
	f1g2 = (atPoint(f1,vars,poss)*atPoint(g2,vars,poss));
	f2g1 = (atPoint(f2,vars,poss)*atPoint(g1,vars,poss));
	
	nextY = poss(2) - ((f1g - fg1)/(f1g2 - f2g1));
end

function value = atPoint(f,vars,pos)
	
	value = eval(subs(subs(f,vars(1),pos(1)),vars(2),pos(2)));

end

