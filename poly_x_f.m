function [ g ] = poly_x_f( f, x, n )
% System of polynomials and polynomial times multiplier function
%   f  Multiplier function
%   x  Function variable
%   n  Number of functions in the system
%   g  System of functions

% Copyright (c) 2017, Juha Sarmavuori

    for i=1:2:n
        g(i) = x^floor(i/2);
    end
    for i=2:2:n
        g(i) = g(i-1) * f;
    end

end

