%HLL with Linearlized Estimation of Wave Speeds
function [Flux, maxSpeed] = HLin(edgeLength,edgeNormal, ULeft, URight)
gamma=1.4;
%[Physical Flux,Density,Velocity normal to edge, Pressure, Sound speed]
%Orientation of the x-axis of Riemann Problem normal to the edge done
%because of putting nV = uL or uR
[FluxL,rhoL,uL,pL,aL] = StateToPFlux(ULeft,edgeNormal,gamma);
[FluxR,rhoR,uR,pR,aR] = StateToPFlux(URight,edgeNormal,gamma);
rhoBar=0.5*(rhoL+rhoR);
aBar=0.5*(aL+aR);
pBar=0.5*(pL+pR);
uBar=0.5*(uL+uR);
pS= pBar - (0.5*(uR-uL)*aBar*rhoBar);
uS= uBar - 0.5*(pR-pL)/(rhoBar*aBar);
rhoSL=rhoL+(((uL-uS)*rhoBar)/aBar);
rhoSR=rhoR+(((uS-uR)*rhoBar)/aBar);
aSL=sqrt(gamma*pS/rhoSL);
aSR=sqrt(gamma*pS/rhoSR);
sL=min([0,uL-aL, uS-aSL]); % making sL zero helps to avoid if else conditions
sR= max([0,uR+aR, uS+aSR]);
Flux=0.5*(FluxL + FluxR)-0.5*(sR+sL)/(sR - sL)*(FluxL - FluxR)+sR*sL/(sR - sL)*(ULeft-URight);
Flux= edgeLength*Flux;
maxSpeed =edgeLength* max(abs(uL) + aL, abs(uR) + aR);