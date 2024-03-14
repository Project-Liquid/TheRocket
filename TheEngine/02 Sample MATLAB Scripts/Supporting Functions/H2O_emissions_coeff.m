function [c,n] = H2O_emissions_coeff(T)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

H2O = [1000 0.165   0.45;
       2000 0.9     0.65;
       3000 2.05    0.61];

c_coeff = polyfit(H2O(:,1),H2O(:,2),2);
c = c_coeff(1)*T.^2+c_coeff(2)*T+c_coeff(3);

n_coeff = polyfit(H2O(:,1),H2O(:,3),2);
n = n_coeff(1)*T.^2+n_coeff(2)*T+n_coeff(3);

end

