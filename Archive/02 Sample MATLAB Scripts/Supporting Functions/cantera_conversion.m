function [output] = cantera_conversion(input)
% Converts a propellant character string compatible with CoolProp into one
% compatible with Cantera. Index of element in CoolProp-compatible input array corresponds with
% index of Cantera-friendly equivalent in output array. Expand for futher
% propellants as you desire

input_array = {'Ethane','Dodecane'};
output_array = {'C2H6','C12H26'};

if ismember(input,input_array) == true % Checks if input has a conversion
    index = find(strcmp(input_array,input));
    output = char(output_array(index));
else % Outputs input if no conversion is found.
    output=char(input);
end

end

