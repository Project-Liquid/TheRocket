function [zero_estimate,value,ea,iteration_number] = bisect_method(f,xl,xu,et)
% bisect_method: bisection root locator that finds roots between xL and xu
% Inputs:
%   f - function to be evaluated
%   xl, xu - lower and upper bounds, respectively
%   et - maximum allowable error (default 0.01% or 0.0001)
% Outputs:
%   zero_estimate – estimated root location
%   value – input function value at the estimated root location
%   error – magnitude of approximate relative error (%)
%   iteration_number - number of iterations the functions goes through

%   Created by: Austin Morse
%   Today's Date: 2/7/21

if nargin < 3, error('At least 3 input arguments required'), end
if nargin == 3
   et = 0.0001;
end

xr=xl;
ea=100;
iteration_number=0;

while ea > et
    iteration_number=iteration_number+1;
    xrprev=xr;
    xr=(xl+xu)/2;
    if f(xl)*f(xr) < 0
        xu=xr;
        ea=abs(xr-xrprev)/xr*100;
    elseif f(xl)*f(xr) > 0
        xl=xr;
        ea=abs(xr-xrprev)/xr*100;
    else
        ea=0;
    end
end

value=f(xr);
zero_estimate=xr;

end

