%*****************************************************************************
%          PRACTICUM NUMERIEKE STROMINLGSLEER 1990
%
% OPGAVE:    1b
% NAAM:      ns_1b.m (volledig)
% DATUM:     14 - 5 - 1990
%*****************************************************************************
% OMSCHRIJVING:
%   Deze module bepaalt een discrete oplossing van
%              dy/dx - k d2y/dx2 = max(2-5x,0)
%   op het interval 0<x<1, met randvoorwaarden y(0)=0 en y(1)=1.
%   Als discretisatie methode wordt een lambda-schema gebruikt.
%*****************************************************************************
% INVOER:   N         aantal roosterpunten
%           k         diffusiecoefficient
%           lambda    coefficient in schema
%                     speciale waarden: 0   centraal
%                                       1/8 QUICK 
%                                       1/2 B3                          
% UITVOER:  ylam(1:N) discrete oplossing
%           A(1:N-1,1:N-1) coefficientenmatrix
%*****************************************************************************
clear upp dia low lowlow A b ylam   % oude arrays verwijderen
%
h=1/N;                              % maaswijdte uniform rooster
%
% coefficienten uitrekenen
%
upp=(0.5/h-lambda/h-k/(h^2))*ones(1,N-1);
dia=(2*k/(h^2)+3*lambda/h)*ones(1,N-1);
low=(-0.5/h-3*lambda/h-k/(h^2))*ones(1,N-2);
lowlow=lambda/h*ones(1,N-3);
%
% eerste vergelijking centraal
%
upp(1)=0.5/h-k/(h^2);
dia(1)=2*k/(h^2);
%
% matrix vullen
%
A=diag(upp(1:N-2),1)+diag(dia)+diag(low,-1)+diag(lowlow,-2);
%
% rechterlid vullen
%
M=round(0.4/h);
b=2-5*h*(1:1:M);
b(M+1:N-1)=zeros(1,N-M-1);
b(N-1)=b(N-1)-upp(N-1);         % rechter randvoorwaarde invullen
%
% oplossing bepalen
%
ylam=inv(A)*b';
ylam(N)=1;
