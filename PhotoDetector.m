classdef PhotoDetector < handle
    properties
        WtoA = 0.62;
        Gain = 30000; %unit: V/A
    end
    
    methods
        function obj = PhotoDetector(WtoA, Gain)
            obj.Gain = Gain;
            obj.WtoA = WtoA;
        end
        
        function SetGain(obj, Gain)
            obj.Gain = Gain;
        end
    end
end