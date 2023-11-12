classdef rocketProperties
    %TANK Tank class with relevant properties
    %Constructor format:
    %
    %   object = rocketProperties()
    %
    %   Object properties must be explicitly defined, i.e.:
    %
    %   obj.property = value;
    %
    %Object Properties:
    %   
    %Stagnation & Chamber Properties
    %pc:          Chamber Pressure       [Pa]
    %Tc:          Chamber Temperature    [K]
    %gamma:       Specific Heat Ratio    [Unitless]
    %R:           Specific Gas Constant  [J/kg-K]
    %         
    %Mass Flow Rates
    %mdot:        Total Mass Flow Rate   [kg/s]
    %mdot_ox:     Ox Mass Flow Rate      [kg/s]
    %mdot_fuel:   Fuel Mass Flow Rate    [kg/s]
    %mdot_BLC:    Coolant Mass Flow Rate [kg/s]
    %         
    %Wall Characteristic Curves
    %L_array:     Contour X-Axis         [m]
    %p_array:     Contour Pressure       [Pa]
    %Tr_array:    Contour Recovery Temp  [K]
    %Tw_array:    Contour Wall Temp      [K]
    %Mach_array:  Axial Mach Array       [Unitless]
    %h_array:     Film Coefficient Array [J/m^2-K]
    %         
    %Rocket Performance & Outputs
    %Thrust:      Thrust                 [N]
    %cstar:       Char. Velocity         [m/s]
    %Cf:          Thrust Coefficient     [Unitless]
    %pe:          Exit Pressure          [Pa]
    %Tmax:        Max temperature        [K]
    %Tthroat:     Throat temperature     [K]
    
    properties
        %Stagnation & Chamber Properties
        pc          %Chamber Pressure       [Pa]
        Tc          %Chamber Temperature    [K]
        gamma       %Specific Heat Ratio    [Unitless]
        R           %Specific Gas Constant  [J/kg-K]
        
        %Mass Flow Rates
        mdot        %Total Mass Flow Rate   [kg/s]
        mdot_ox     %Ox Mass Flow Rate      [kg/s]
        mdot_fuel   %Fuel Mass Flow Rate    [kg/s]
        mdot_BLC    %Coolant Mass Flow Rate [kg/s]
        
        %Wall Characteristic Curves
        L_array     %Contour X-Axis         [m]
        p_array     %Contour Pressure       [Pa]
        Tr_array    %Contour Recovery Temp  [K]
        Tw_array    %Contour Wall Temp      [K]
        Mach_array  %Axial Mach Array       [Unitless]
        h_array     %Film Coefficient Array [J/m^2-K]
        
        %Rocket Performance & Outputs
        thrust      %Thrust                 [N]
        ISP         %Specific Impulse       [s]
        cstar       %Char. Velocity         [m/s]
        Cf          %Thrust Coefficient     [Unitless]
        pe          %Exit Pressure          [Pa]
        Tmax        %Max temperature        [K]
        Tthroat     %Throat temperature     [K]
    end
    
    methods
        function obj = rocketProperties()
            %ROCKETPROPERTIES Construct an instance of this class.
            %Properties must be manually defined i.e.:
            %
            %obj.property = value;
        end
    end
end

