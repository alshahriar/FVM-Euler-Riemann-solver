function [PFlux,rho,nV,p,a] = StateToPFlux(U,edgeN,g)
rho = U(1); u = U(2)/rho; v = U(3)/rho; E=U(4);
pV=sqrt(U(2)^2 + U(3)^2)/rho; % Physical Velocity
nV=edgeN(1)*u+edgeN(2)*v; %Normal Velocity
%tV = sqrt(pV^2 - nV^2); % Tangential Velocity
p = (g-1)*(E - 0.5*rho*pV^2);
a = sqrt(g*p/(rho));
PFlux = [rho*nV; ...
rho*nV*u + p*edgeN(1); ...
rho*nV*v + p*edgeN(2); ...
(E + p)*nV];