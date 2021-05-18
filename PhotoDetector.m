classdef PhotoDetector < handle
    properties
        WtoA = 0.5;
        Gain = 30;
    end
    
    methods
        function obj = PhotoDetector(Gain)
            obj.Gain = Gain;
        end
        
        function SetGain(obj, Gain)
            obj.Gain = Gain;
        end
    end
end