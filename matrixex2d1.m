% Copyright (c) 2017, Juha Sarmavuori

x=sym('x');
y=sym('y');
f = [1 x+y x*y (x+y)^2 (x*y)^2 (x+y)^3 (x*y)^3 (x+y)^4 (x*y)^4 (x+y)^5 (x*y)^5 (x+y)^6 (x*y)^6 (x+y)^7 (x*y)^7 (x+y)^8 (x*y)^8  (x+y)^9 (x*y)^9];

g1 = x*y;
g2 = x+y;
g = [g1 g2];

M = mult_sym(g, f, 1, [x y], [0 0], [1 1])
M = double(M);

clear app;
for i=1:length(f)
    M3 = expm(M(1:i,1:i,1))*logm(eye(i)+M(1:i,1:i,2));
    app(i) = M3(1,1);
end

[exact, err] = quad2d(@(x,y) exp(x.*y).*log(1+x+y),0,1,0,1)
replace(sprintf('%d & $%s$ & %.15E & %.15E \\\\\n', [1:length(f);f;app;exact-app]),'*',' \, ')
