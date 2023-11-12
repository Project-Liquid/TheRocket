function Output = CoolProp(OutputVar,InputVar1,InputVal1,InputVar2,InputVal2,Fluid)
    Output = py.CoolProp.CoolProp.PropsSI(OutputVar,InputVar1,InputVal1,InputVar2,InputVal2,Fluid);
end
% This just shortens the function to call CoolProp.
% You must have CoolProp installed to use this.

% Input/Outputs
%   'P' : Pressure (Pa)
%   'T' : Temperature (K)
%   'HMOLAR' : Molar Specific Enthalpy (J/mol)
%   'HMASS" : Mass Specific Enthalpy (J/Kg)
%   'UMOLAR' : Molar specific internal energy (J/mol)
%   'U' : Mass Specific Internal Energy (J/Kg)
%   'DMOLAR' : Molar Density (mol/m^3)
%   'D' : Mass density (kg/m^3)
%   'SMOLAR' : Molar specific entropy (J/mol/K)
%   'S' : Mass specific entropy (J/kg/K)
%   'CVMOLAR' : Molar specific constant volume specific heat (J/mol/K)
%   'CVMASS' : Mass specific constant volume specific heat (J/Kg/K)
%   'Cpmass' : Mass specific constant pressure specific heat (J/Kg/K)
%   'PCRIT' : Pressure at the critical point (Pa)
%   'PRANDTL' : Prandtl number
%   'V' : Viscosity (Pa s)
%   'Z' : Compressibility factor