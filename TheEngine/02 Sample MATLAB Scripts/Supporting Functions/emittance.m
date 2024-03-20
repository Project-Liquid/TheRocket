function [em] = emittance(exhaust,Aw,D)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

T_g = temperature(exhaust);
p_g = pressure(exhaust)/101325;

H2O_mol_frac = massFraction(exhaust,'H2O')*meanMolecularWeight(exhaust)/(CoolProp('M','T',T_g,'P',p_g,'H2O')*1000);
CO2_mol_frac = massFraction(exhaust,'CO2')*meanMolecularWeight(exhaust)/(CoolProp('M','T',T_g,'P',p_g,'CO2')*1000);

H2O_p = p_g * H2O_mol_frac;
CO2_p = p_g * CO2_mol_frac;

L_eff = 0.95*D*Aw^-0.85;

H2O_optD = H2O_p * L_eff;
CO2_optD = CO2_p * L_eff;

optD = H2O_optD+CO2_optD;

[c_H2O,n_H2O] = H2O_emissions_coeff(T_g);

[c_CO2,n_CO2] = CO2_emissions_coeff(T_g);

C1 = 0.26 + 0.74*exp(-2.5*H2O_optD);
C2 = 0.75 + 0.31*exp(-10*H2O_optD);
K_H2O = 1 + C1*(1-exp((1-p_g*(1+H2O_mol_frac))/C2));

em_H2O = K_H2O*0.825*(1+(optD/c_H2O)^-n_H2O)^(-1/n_H2O);

m = 100*CO2_optD;
K_CO2 = 10^(0.036*CO2_optD^(-0.433)*(1+(2*log10(p_g))^-m)^(-1/m));

em_CO2 = K_CO2*0.231*(1+(optD/c_CO2)^-n_CO2)^(-1/n_CO2);

n = 5.5*(1+(1.09*optD)^-3.33)^(-1/3.33);
K = 1 - abs((2*H2O_mol_frac/(CO2_mol_frac+H2O_mol_frac))-1)^n;
em_correction = 0.0551*K*(1-exp(-4*optD))*(1-exp(-12.5*optD));


em = em_H2O+em_CO2-em_correction;
end

