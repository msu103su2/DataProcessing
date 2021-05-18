classdef SplitPair < handle
    properties
        d;
        offset;
    end
    methods
        function obj = SplitPair(gapsize, offset)
            obj.d = gapsize;
            obj.offset = offset;            
        end
        
        function result = Fweight(obj, x, y)
            result = heaviside(x - obj.offset - obj.d/2) - heaviside(-(x - obj.offset + obj.d/2));
        end
        
        function result = F_onlyx(obj, HGBeam)
            scale = 3;
            waypoints = -obj.offset:0.1*obj.offset:obj.offset;
            func = @(x) HGBeam.u(x, x, 0) .* obj.Fweight(x, 0);
            result = integral(func, -scale*HGBeam.w0x, scale*HGBeam.w0x,'RelTol',1e-3,'Waypoints',waypoints);
        end
    end
end