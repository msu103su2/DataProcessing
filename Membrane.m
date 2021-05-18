classdef Membrane < handle
    properties
        a;
        b;
        h;
        rho;
        s;
        fm;
        gamma;
        modeIndex;
    end
    
    methods
        function obj = Membrane(a,b,h,rho)
            obj.a = a;
            obj.b = b;
            obj.h = h;
            obj.rho = rho;
        end
        
        function freqs = modes(obj, c)
            Is = [];
            Js = [];
            freqs = [];
            for i = 1:100
                for j = 1:100
                    freqs = [freqs, c * sqrt((i*pi/obj.a)^2 + (j*pi/obj.b)^2)];
                    Is = [Is, i];
                    Js = [Js, j];
                end
            end
            freqs = [freqs; Is; Js];
            freqs = freqs';
            freqs = sortrows(freqs);
        end
        
        function SetStress(obj, stress)
            obj.s = stress;
        end
        
        function SetMode(obj, fm, gamma, modeIndex)
            obj.fm = fm;
            obj.gamma = gamma;
            obj.modeIndex = modeIndex;
        end
        
        function chi = response(obj, f)
            omegaM = 2*pi*obj.fm;
            omega = 2*pi*f;
            meff = obj.rho*obj.a*obj.b*obj.h/4;
            chi = 1/meff*((omegaM^2 - omega.^2)+1i*obj.gamma*omega).^(-1);
        end
        
        function re = kmx(obj)
            re = obj.modeIndex*pi/obj.a;
        end
        
        function m = meff(obj)
            m = obj.rho*obj.a*obj.b*obj.h/4;
        end
    end
end