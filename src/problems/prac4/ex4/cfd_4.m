%***************************************************************************
%           JMBC PRACTICUM NUMERIEKE STROMINGSLEER 1998
%
% OPGAVE :   4
% NAAM MODULE: cfd_4.m
% DATUM : January 1998 
%***************************************************************************
% OMSCHRIJVING:
%   Deze module lost een warmtetransportvergelijking op
%            dT/dt + u dT/dx = k d2T/dx2
%   voor 0<x<1 en 0<T<4. Randvoorwaarden zijn T(0,t)=0 en dT/dx(1,t)=0.
%   Beginvoorwaarde T(x,0)=x.
%   Discretisatie methode is gegen. Crank-Nicolson met parameter omega.
%
%**************************************************************************
% INVOER:   NX       aantal roosterpunten in x-richting
%           dt       tijdstap
%           omega    gewichtsfactor in tijdsintegratie
%           u        convectieve snelheid
%           k        diffusiecoefficient 
% UITVOER:  tijd(i)  tijd t=(i-1)*dt, i=1:NTMAX+1
%           temp(i)  T(1,tijd(i))      temperatuur op x=1
%           T(i,n)   T(x(i),tijd(n))   volledige oplossing
%***************************************************************************
%
clear                      % alle variabelen verwijderen 
%
% begin invoer
%
%u=input('input velocity u:   ');
u=0;
%k=input('input diffusion k:    ');
k=1;
NX=input('input grid points N:   ');
omega=input('input parameter omega:   ');
dt=input('input time step dt:   ');
tmax=4;                       % lengte tijdsinterval
%
% einde invoer 
%
NTMAX=round(tmax/dt);         % aantal tijdstappen
hx=1/NX;                      % maaswijdte
d=2*k*dt/(hx*hx); 
eta = u*dt/hx;                % CFL getal
%
% berekenen matrix coefficienten
%
for i=1:NX-1,
   di(i)=1+omega*d;
   lo(i)=omega*0.5*(-eta-d);
   up(i)=omega*0.5*(eta-d);
end
%
% op x=1 een Neumann randvoorwaarde
%
di(NX)=1+omega*d;
lo(NX)=-omega*d;
up(NX)=0;
%
% vullen matrix 
%
A=diag(up(1:NX-1),1)+diag(di)+diag(lo(2:NX),-1);
% eenmalig A inverteren
AINV=inv(A);
%
% beginvoorwaarde invullen
%
x=[1:NX]*hx;
tnew=x;
%
T=tnew;
tijd(1)=0;                    % begintijdstip
temp(1)=1;                    % beginvoorwaarde
%
% uitvoeren tijdstappen -------------------------------------------
%
for NT=1:NTMAX,                % start tijdstappen
   told=tnew;
%
% berekenen rechterlid
%
   rl(1)=(1-omega)*0.5*(-2*d*told(1)+(d-eta)*told(2)) + told(1);
%
   for i=2:NX-1,
      rl(i)=(1-omega)*0.5*((d+eta)*told(i-1)-2*d*told(i)+(d-eta)*told(i+1)) + told(i);
   end
%
% op x=1 een Neumann randvoorwaarde
%
   rl(NX)=(1-omega)*d*(told(NX-1)-told(NX))+told(NX);
%
% oplossen
%
   tnew=AINV*rl';
%
%antwoord opslaan
%
   tijd(NT+1)=NT*dt;            % nummering is een versprongen
   temp(NT+1)=tnew(NX);        % nummering is een versprongen

   %print output to screen
   if (dt<=0.1)
    for i=1:40
     if ((tijd(NT+1)<=i/10+dt/2) & (tijd(NT+1)>=i/10-dt/2))
      if (abs(temp(NT+1)) <=1)
       fprintf('t=%1.2f   Temp(1,t)= %1.4f \n',tijd(NT+1),temp(NT+1))
      else
       fprintf('t=%1.2f   Temp(1,t)= %4.2e \n',tijd(NT+1),temp(NT+1))    
      end
     end
    end
   else
    if (abs(temp(NT+1)) <=1)
       fprintf('t=%1.2f   Temp(1,t)= %1.4f \n',tijd(NT+1),temp(NT+1))
    else
       fprintf('t=%1.2f   Temp(1,t)= %4.2e \n',tijd(NT+1),temp(NT+1))    
    end
   end
   
   T=[T;[tnew(1:NX)']];
end                           % einde tijdstappen   
%
%-----------------------------------------------------------------
%                             
clf                           % begin plaatjes
subplot(121)
plot(tijd,temp)
set(gca,'XLim',[0 tmax]);
title('solution at x=1')
xlabel('t')
ylabel('Temp')
tekst=sprintf('dt=%6.5f',dt);
text(0.1,0.9,tekst,'sc')
tekst=sprintf('omega=%4.2f',omega);
text(0.5,0.9,tekst,'sc')
%
subplot(122)
mesh(tijd,[0 x],[zeros(size(tijd,2),1) T]')
set(gca,'XLim',[0 tmax]);
title('Temp(x,t)')
xlabel('t')
ylabel('x')
zlabel('Temp')
view(30,30)







