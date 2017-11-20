% Copyright (c) 2017, Juha Sarmavuori

% Generalized Gaussian quadrature nodes and weights.
% Taken from table 2 of:
% Generalized Gaussian Quadrature Rules for Systems of Arbitrary Functions
% by J. Ma, V. Rokhlin and S. Wandzura
% in SIAM Journal on Numerical Analysis, Vol. 33, No. 3 (Jun., 1996), pp. 971-996
% or
% table 6 of:
% Generalized Gaussian Quadrature Rules for Systems of Arbitrary Functions
% by J. Ma, V. Rokhlin and S. Wandzura
% in Research Report YALEU/DCS/RR-990 October, 1993
nw5 = [
0.831174531456776E-02 0.279778782123048E-01
0.863966362795308E-01 0.142398935990482E+00
0.305943516943443E+00 0.291943668689807E+00
0.635656558720652E+00 0.339447240627363E+00
0.920004207024857E+00 0.198232276480043E+00];

nw20 =[
0.565910606963981E-04 0.196275471368437E-03
0.680854395380959E-03 0.125263851254803E-02
0.306270110956997E-02 0.380561564887333E-02
0.894917200513315E-02 0.832217213031223E-02
0.204442383726108E-01 0.150414196754775E-01
0.397611592475325E-01 0.239388575134562E-01
0.689512869857525E-01 0.347155654274889E-01
0.109635234283660E+00 0.468137478429227E-01
0.162762389006638E+00 0.594574979240510E-01
0.228422459407150E+00 0.717153419416687E-01
0.305728258549570E+00 0.825791227715473E-01
0.392782689780531E+00 0.910522631915573E-01
0.486735430814252E+00 0.962395518442346E-01
0.583926798874029E+00 0.974303994336976E-01
0.680108431067843E+00 0.941680397096348E-01
0.770723438587922E+00 0.862983560861100E-01
0.851223225892569E+00 0.739938065555806E-01
0.917394745239490E+00 0.577502214287376E-01
0.965671257181204E+00 0.383577114718547E-01
0.993407695194012E+00 0.168713954188783E-01];

x = sym('x');
n = 5;

f = poly_x_f(x^(1/3), x, n);

g = [x x^sym(1/3)];

Mx = mult_sym(g, f, 1, x, 0, 1);

[v1,d1] = eig(double(Mx(:,:,1)));
nodes1 = diag(d1);
weights1 = abs(v1(1,:)).^2';

[v2,d2] = eig(double(Mx(:,:,2)));
nodes2 = diag(d2).^3;
weights2 = abs(v2(1,:)).^2';

ssp1=sortrows([nodes1 weights1], 1);
ssp2=sortrows([nodes2 weights2], 1);

scatter(ssp1(:,1),ssp1(:,2),'k+');
hold on;
scatter(ssp2(:,1), ssp2(:,2),'kx');
scatter(nw5(:,1),nw5(:,2),'ko')
legend('Matrix method g(x)=x', 'Matrix Method g(x)=x^{1/3}', 'Generalized Gaussian', 'Location', 'south');
xlabel('node x_i');
ylabel('weight w_i');
%matlab2tikz('ex1_nodes.tex', 'standalone', true);
hold off

ys = sym('y');
powers = [0:0.01:6.5];
exact = 0*powers;
gg = 0*powers;
m1 = 0*powers;
m2 = 0*powers;
pint = int(x^sym(ys), x, 0, 1);

for i=1:length(powers)
    exact(i) = double(subs(pint,ys,powers(i)));
    gg(i) = sum(nw5(:,2) .* nw5(:,1) .^ powers(i));
    m1(i) = sum(ssp1(:,2) .* ssp1(:,1) .^ powers(i));
    m2(i) = sum(ssp2(:,2) .* ssp2(:,1) .^ powers(i));
end

plot(powers, m1 ./ exact - 1, 'k:', powers, m2 ./ exact - 1, 'k--', powers, gg ./ exact - 1, 'k-');
legend('Matrix method g(x)=x', 'Matrix method g(x)=x^{1/3}', 'Generalized Gaussian', 'Location', 'northEast');
ylabel('Error: \epsilon');
xlabel('Power: y');
%matlab2tikz('ex1_error5.tex', 'standalone', true,'extraAxisOptions','scale=\figurescale');

n = 20;

f = poly_x_f(x^(1/3), x, n);

g = [x x^sym(1/3)];

Mx = mult_sym(g, f, 1, x, 0, 1);

[v1,d1] = eig(double(Mx(:,:,1)));
nodes1 = diag(d1);
weights1 = abs(v1(1,:)).^2';

[v2,d2] = eig(double(Mx(:,:,2)));
nodes2 = diag(d2).^3;
weights2 = abs(v2(1,:)).^2';

ssp1=sortrows([nodes1 weights1], 1);
ssp2=sortrows([nodes2 weights2], 1);

ys = sym('y');
%powers = [0:0.01:30];
exact = 0*powers;
gg = 0*powers;
m1 = 0*powers;
m2 = 0*powers;
pint = int(x^sym(ys), x, 0, 1);

for i=1:length(powers)
    exact(i) = double(subs(pint,ys,powers(i)));
    gg(i) = sum(nw20(:,2) .* nw20(:,1) .^ powers(i));
    m1(i) = sum(ssp1(:,2) .* ssp1(:,1) .^ powers(i));
    m2(i) = sum(ssp2(:,2) .* ssp2(:,1) .^ powers(i));
end

plot(powers, m1 ./ exact - 1, 'k:', powers, m2 ./ exact - 1, 'k--', powers, gg ./ exact - 1, 'k-');
legend('Matrix method g(x)=x', 'Matrix method g(x)=x^{1/3}', 'Generalized Gaussian', 'Location', 'northEast');
ylabel('Error: \epsilon');
xlabel('Power: y');
%matlab2tikz('ex1_error20.tex', 'standalone', true, 'extraAxisOptions', 'scale=\figurescale');
