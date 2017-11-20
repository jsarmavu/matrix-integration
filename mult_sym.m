function [ M ] = mult_sym( g, f, w, x, a, b )
% Matrix approximation of multiplication operator
%   g  Functions for the multiplication operator
%   f  Linearly independent functions
%   w  Weight function
%   x  Integration variable
%   a  Lower limit of integration interval
%   b  Upper limit of integration interval
%   M  Matrix approximation of multiplication with function g

% Copyright (c) 2017, Juha Sarmavuori

    n = length(f);

    for m=1:length(g)
        for i=1:n
            for j=i:n
                Mi = f(i)' * w * g(m) * f(j);
                for k=1:length(x)
                    Mi = int(Mi, x(k), a(k), b(k));
                end
                M(i,j,m) = Mi;
            end
        end
    end
    for m=1:length(g)
        M(:,:,m)  = M(:,:,m) + M(:,:,m)' - diag(diag(M(:,:,m)));
    end

    for i=1:n
        for j=i:n
            Gi = f(i)' * w * f(j);
            for k=1:length(x)
                Gi = int(Gi, x(k), a(k), b(k));
            end
            ip(i,j) = Gi;
        end
    end
    ip = ip + ip' - diag(diag(ip));

    R = chol(ip);

    for m=1:length(g)
        M(:,:,m) = R' \ M(:,:,m) / R;
    end
end

