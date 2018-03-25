
classdef Function
    properties
        vars
        func
    end
    
    methods
        
        function obj = Function(f, variables)
            obj.func = f;
            obj.vars = variables;
        end
        
        function pd = parDeriv(Func, pdNum)
            f = diff(Func.func, Func.vars(pdNum));
            pd = Function(f,Func.vars);
        end
        
        function value = valueOf(Func, point)
            sub = subs(Func.func, Func.vars(1), point(1));
            for i = 2:length(Func.vars)
                sub = subs(sub, Func.vars(i), point(i));
            end
            value = eval(sub);
        end
        
    end
end

