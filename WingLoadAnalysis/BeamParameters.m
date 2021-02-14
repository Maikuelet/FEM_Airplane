function [A, Iy, Iz, J] = BeamParameters(mat,dat,Tmat, e)
% Computes Beam section area, inertia and torsional constant
% for an introduced element

m = Tmat(e);

%% Section parameters
a=mat(6,m);
b=mat(7,m);
t=mat(8,m);
h=dat(1,e);

%% Parameter Computation
A   = (h-t)*a+2*b*t;
Iy  = (1/12)*a*(h-t)^3+(1/6)*b*t^3+b*t*(h^2/2);
Iz  = (1/12)*(h-t)*a^3+(1/6)*t*b^3;
J   = (1/3)*(2*b*t^3+h*a^3);

end