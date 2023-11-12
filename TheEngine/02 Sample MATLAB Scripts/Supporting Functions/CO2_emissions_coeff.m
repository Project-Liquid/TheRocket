function [c,n] = CO2_emissions_coeff(T)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

CO2 = [1000 0.05   0.6;
       1500 0.075  0.6;
       2000 0.15   0.6];

c_coeff = polyfit(CO2(:,1),CO2(:,2),2);
c = c_coeff(1)*T.^2+c_coeff(2)*T+c_coeff(3);

n_coeff = polyfit(CO2(:,1),CO2(:,3),2);
n = n_coeff(1)*T.^2+n_coeff(2)*T+n_coeff(3);

end

