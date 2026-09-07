%*****************************************************************************
%     PRACTICUM NUMERIEKE STROMINGSLEER 
%
% OPGAVE: 1a
% NAAM: cfd_1.m
% DATUM: september 2001
%******************************************************************************
% OMSCHRIJVING:
%    Deze module plot 3 upwind oplossingen en  1 exacte oplossing van 
%        du/dx - k d2u/dx2 = max(2 - 5x,0)
%    op het interval 0 < x < 1, met randvoorwaarden y(0)=0 en y(1)=1
%******************************************************************************
% INVOER:   N  aantal roosterpunten
%           k  diffusiecoefficient 
% UITVOER:  plot met 3 upwind en 1 exacte oplossing
%******************************************************************************

clear x yupw yup1 yup2 yup3 yex 
clf

N=input('give value of N:  ');
x=1/N*[1:1:N];                          % definitie uniform rooster
methode='U';

k=input('give first value for k:  ' );
cfd_upw
yup1=yupw;
plot(x,yup1,'r+')
hold on
tekst=sprintf(' k=%6.4f',k);
plot([0.15 0.18],[0.6 0.6],'r+')
text(0.2,0.6,tekst)
title('upwind computations')
xlabel('x')
axis([0 1 0 1]);

k=input('give second value for k:  ');
cfd_upw
yup2=yupw;
plot(x,yup2,'bo')
tekst=sprintf(' k=%6.4f',k);
plot([0.15 0.19],[0.7 0.7],'bo')
text(0.2,0.7,tekst)

k=input('give third value for k:  ');
cfd_upw
yup3=yupw;
plot (x,yup3,'m-.')
tekst=sprintf(' k=%6.4f',k);
plot([0.15 0.19],[0.8 0.8],'m-.')
text(0.2,0.8,tekst)

k=input('give value for k of exact solution:   ');
for i=1:N
   yex(i)=ns_exact(x(i),k);             % berekening exacte oplossing
end
plot(x,yex,'g-')
tekst=sprintf(' k=%6.4f',k);
plot([0.15 0.19],[0.9 0.9],'g-')
text(0.2,0.9,[tekst ' exact'])
