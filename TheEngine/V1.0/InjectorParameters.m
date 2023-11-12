classdef InjectorParameters < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Cd;
        elementDiam;
        elementCount;
        elementArea;
        cumulativeArea;
    end
    
    methods
        function obj = InjectorParameters(Cd, elementDiam, elementCount)
            obj.Cd = Cd;
            obj.elementDiam = elementDiam;
            obj.elementCount = elementCount;
            obj.elementArea = (pi/4) * (elementDiam^2);
            obj.cumulativeArea = obj.elementArea * elementCount;
        end
    end
end

