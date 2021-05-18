classdef Laser < handle
    properties
        lambda;
        Pdc;
        k;
        f;
        omega;
        Pac;
        deltaf;
    end
    
    methods
        function obj = Laser(lambda)
            temp = consts;
            obj.lambda = lambda;
            obj.k = 2*pi/lambda;
            obj.f = temp.c/lambda;
            obj.omega = 2*pi*obj.f;
        end
        
        function SetPdc(obj, Pdc)
            obj.Pdc = Pdc;
        end
        
        function SetAC(obj, Pac, deltaf)
            obj.Pac = Pac;
            obj.deltaf = deltaf;
        end
    end
end