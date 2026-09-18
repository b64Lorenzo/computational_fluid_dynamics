%*****************************************************************************
%          PRACTICUM NUMERIEKE STROMINLGSLEER 1993
%
% OPGAVE:    1c
% NAAM:      ns_opg1c.m
% DATUM:     maart 1993
%*****************************************************************************
% OMSCHRIJVING:
%   Deze module bepaalt een discrete oplossing van
%              dy/dx - k d2y/dx2 = max(2-5x,0)
%   op het interval 0<x<1, met randvoorwaarden y(0)=0 en y(1)=1.
%   Er wordt centraal of upwind gediscretiseerd op een gerekt rooster.
%*****************************************************************************
% INVOER:   N         aantal roosterpunten
%           k         diffusiecoefficient
%           methode   methode A, B of U
%
% UITVOER:  xrek(1:N)    roosterpunten
%           yexrek(1:N)  exacte oplossing 
%           yrek(1:N)    discrete oplossing
%*****************************************************************************
clear xrek yrek upp dia low A b  % oude arrays weggooien 
%
% berekening van gerekte rooster; N/2 mazen links en rechts van x = 1-Lk
% goede keuze L=5
L=5;
%
hcrs=2.*(1.-L*k)/N;           % grove mazen links
hfin=2*L*k/N;                 % fijne mazen rechts
for i=1:N/2
   xrek(i)=i*hcrs;            % roosterpunten links
end
for i=N/2+1:N 
   xrek(i)=1.-(N-i)*hfin;     % roosterpunten rechts
end
% 
% exacte oplossing op gerekt rooster 
%
for i=1:N
yexrek(i)=ns_exact(xrek(i),k);
end
%
if methode=='U'
%
% coefficienten uitrekenen upwind
%
   upp(1:N/2-1)=-k/(hcrs^2)*ones(1,N/2-1);
   dia(1:N/2-1)=(1.0/hcrs+2.*k/(hcrs^2))*ones(1,N/2-1);
   low(1:N/2-1)=(-1.0/hcrs-k/(hcrs^2))*ones(1,N/2-1);
%		  
   upp(N/2)=-2.*k/(hfin*(hfin+hcrs));
   dia(N/2)= 1.0/hcrs+2.*k/(hfin*hcrs);
   low(N/2)=-1.0/hcrs-2.*k/(hcrs*(hcrs+hfin));
%
   upp(N/2+1:N-1)=-k/(hfin^2)*ones(1,N/2-1);
   dia(N/2+1:N-1)=(1.0/hfin+2.*k/(hfin^2))*ones(1,N/2-1);
   low(N/2+1:N-1)=-(1.0/hfin+k/(hfin^2))*ones(1,N/2-1);
else
%
% coefficienten uitrekenen centraal
%
   upp(1:N/2-1)=(0.5/hcrs-k/(hcrs^2))*ones(1,N/2-1);
   dia(1:N/2-1)=2.*k/(hcrs^2)*ones(1,N/2-1);
   low(1:N/2-1)=-0.5/hcrs-k/(hcrs^2)*ones(1,N/2-1);
%
   i=N/2;
   if methode=='A' 
      upp(i)=1./(hcrs+hfin)-2.*k/(hfin*(hfin+hcrs));
      dia(i)=2.*k/(hfin*hcrs);
      low(i)=-1./(hcrs+hfin)-2.*k/(hcrs*(hcrs+hfin));
   end
   if methode=='B' 
      upp(i)=(hcrs-2.*k)/(hfin*(hfin+hcrs));
      dia(i)=(hfin-hcrs+2.*k)/(hfin*hcrs);
      low(i)=(-hfin-2.*k)/(hcrs*(hcrs+hfin));
   end
%
   upp(N/2+1:N-1)=(0.5/hfin-k/(hfin^2))*ones(1,N/2-1);
   dia(N/2+1:N-1)=2.*k/(hfin^2)*ones(1,N/2-1);
   low(N/2+1:N-1)=(-0.5/hfin-k/(hfin^2))*ones(1,N/2-1);
end
%
% matrix vullen
%
A=diag(upp(1:N-2),1)+diag(dia(1:N-1),0)+diag(low(2:N-1),-1);
%
% rechterlid vullen
%
M=round(0.4/hcrs);
b=2.-5.*hcrs*(1:1:M);
b(M+1:N-1)=zeros(1,N-M-1);
b(N-1)=b(N-1)-upp(N-1);       % randvoorwaarde rechts invullen
%
% oplossing bepalen
%
yrek=A\b';
yrek(N)=1;


