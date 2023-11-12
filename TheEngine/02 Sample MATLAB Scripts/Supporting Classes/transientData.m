classdef transientData
    %Transient data handling class with relevant properties. Same format as
    %
    %Constructor format:
    %
    %   object = transientData()
    %
    %   Object properties must be explicitly defined, i.e.:
    %
    %   obj.property = value;
    %
    %Object Properties:
    %   
    %Tank Conditions
    %T0Ox        Ox tank temperature    [K]
    %p0Ox        Ox tank pressure       [Pa]
    %QOx         Ox drain percentage    [%]
    %mOx         Ox mass                [kg]
    %         
    %T0Fuel      Fuel tank temperature  [K]
    %p0Fuel      Fuel tank pressure     [Pa]
    %QFuel       Fuel drain percentage  [%]
    %mFuel       Fuel mass              [kg]
    %         
    %T0BLC       BLC tank temperature   [K]
    %p0BLC       BLC tank pressure      [Pa]
    %QBLC        BLC drain percentage   [%]
    %mBLC        BLC mass               [kg]
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
    %thrust:      Thrust                 [N]
    %cstar:       Char. Velocity         [m/s]
    %Cf:          Thrust Coefficient     [Unitless]
    %pe:          Exit Pressure          [Pa]
    %Tmax:        Max temperature        [K]
    %Tthroat:     Throat temperature     [K]
    
    properties
        %Transient Tank Conditions
        T0Ox        %Ox tank temperature    [K]
        p0Ox        %Ox tank pressure       [Pa]
        QOx         %Ox drain percentage    [%]
        mOx         %Ox mass                [kg]
        
        T0Fuel      %Fuel tank temperature  [K]
        p0Fuel      %Fuel tank pressure     [Pa]
        QFuel       %Fuel drain percentage  [%]
        mFuel       %Fuel mass              [kg]
        
        T0BLC       %BLC tank temperature   [K]
        p0BLC       %BLC tank pressure      [Pa]
        QBLC        %BLC drain percentage   [%]
        mBLC        %BLC mass               [kg]
        
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
        cstar       %Char. Velocity         [m/s]
        Cf          %Thrust Coefficient     [Unitless]
        pe          %Exit Pressure          [Pa]
        Tmax        %Max temperature        [K]
        Tthroat     %Throat temperature     [K]
    end
    
    methods
        function obj = transientData()
            %TRANSIENTDATA Construct an instance of this class.
            %Properties must be manually defined i.e.:
            %
            %obj.property = value;
        end
    end
end

