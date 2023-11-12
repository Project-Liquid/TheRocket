classdef EngineParameters < handle
    %All parameters of the engine
    
    properties
        diamThroat;
        areaThroat;
        Pc;
        m_dot;
        m_dotox;
        m_dotfuel;
        Tc;
        gamma;
        R;
        expRatio;
        charV;
        Cf;
        thrust;
        Isp;
    end
    
    methods
        function obj = EngineParameters(diamThroat, expRatio, i)
            
            %Preallocation
            preallocation = 10000;
            
            obj.Pc = zeros(1,preallocation);
            obj.m_dot = zeros(1,preallocation);
            obj.m_dotox = zeros(1,preallocation);
            obj.m_dotfuel = zeros(1,preallocation);
            obj.Tc = zeros(1,preallocation);
            obj.gamma = zeros(1,preallocation);
            obj.R = zeros(1,preallocation);
            obj.charV = zeros(1,preallocation);
            obj.Cf = zeros(1,preallocation);
            obj.thrust = zeros(1,preallocation);
            obj.Isp = zeros(1,preallocation);
            
            %Set Static Parameters
            obj.diamThroat = diamThroat;
            obj.areaThroat = (pi/4) * (diamThroat^2);
            obj.expRatio = expRatio;
        end
    end
end

