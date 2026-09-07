%*****************************************************************************
%     PRACTICUM NUMERIEKE STROMINGSLEER 
%
% OPGAVE: 1a
% NAAM: ns_1a.m
% DATUM: mei 1993
%******************************************************************************
% OMSCHRIJVING:
%    Deze module berekent de upwind en de smart-upwind oplossing van 
%        du/dx - k d2u/dx2 = max(2 - 5x,0)
%    op het interval 0 < x < 1, met randvoorwaarden y(0)=0 en y(1)=1
%******************************************************************************
% INVOER:   N  aantal roosterpunten
%           k  diffusiecoefficient
%           methode 'U' of 'S'
% UITVOER:  yupw(1:N)  discrete oplossing
%           A(1:N-1,1:N-1) coefficientenmatrix
%******************************************************************************
clear upp dia low A b yupw     % oude arrays verwijderen
%
h=1/N;                         % maaswijdte uniform rooster
%
% coefficienten berekenen
%
if methode=='U'
   upp=-k/(h^2)*ones(1,N-1);               % bovendiagonaal
   dia=(2*k/(h^2)+1/h)*ones(1,N-1);        % diagonaal
   low=(-k/(h^2)-1/h)*ones(1,N-1);         % onderdiagonaal
elseif methode=='S'
   keff=0.5*h/tanh(h/(2*k));
   upp=(0.5/h-keff/(h^2))*ones(1,N-1);
   dia=2*keff/(h^2)*ones(1,N-1);
   low=(-0.5/h-keff/(h^2))*ones(1,N-1);
end
%
% matrix vullen
%
A=diag(upp(1:N-2),1)+diag(dia(1:N-1),0)+diag(low(2:N-1),-1);  
%
% rechterlid berekenen
%
M=round(0.4/h);
b=2-5*h*(1:1:M);
b(M+1:N-1)=zeros(1,N-M-1);
b(N-1)=b(N-1)-upp(N-1);         % rechterrandvoorwaarde invullen
%
% oplossing bepalen
%
yupw=A\b';
yupw(N)=1;

