%******************************************************************************
%               PRAKTIKUM NUMERIEKE STROMINGSLEER 1993
%
% OPGAVE:       1
% NAAM MODULE:  ns_exact.m
% DATUM:        maart 1993
%******************************************************************************
% OMSCHRIJVING:
%    Deze module bepaalt de exacte oplossing van 
%           dy/dx - k d2y/dx2 = max(2 - 5x,0)
%    op het interval 0 < x < 1, met randvoorwaarden y(0)=0, y(1)=1.
%    De oplossing wordt berekend in het opgegeven roosterpunt x.
%******************************************************************************
% INVOER:            x   coordinaat roosterpunt
%                    k   diffusiecoefficient
% UITVOER:      yexact   exacte oplossing in roosterpunt
%******************************************************************************
%
function yexact=ns_exact(x,k)
%
pe=1/k;
% aktie om underflow van exp(-pe) te voorkomen
if pe > 100,
   empe=0;
else
   empe=exp(-pe);
end
%
% algemeen rechterlid is rl*(x1-x)
rl=5;
x1=0.4;
%
% op x < x1 is oplossing a1*x^2 + b1*x + c1 + d1*exp(pe*(x-x1))
% op x > x1 is oplossing c2 + d2*exp(pe*(x-1))
%
a1=-0.5*rl;
b1=(0.4*rl+2*k*a1);
f1=a1*x1*x1+b1*x1-(2*x1*a1+b1)/pe;
g1=(2*a1*x1+b1)/pe;
d2=(1-g1*exp(-pe*x1)-f1)/(1-empe);
c2=1-d2;
c1=c2-f1;
d1=d2*exp(pe*(x1-1))-g1;
%
% oplossing kan nu uitgerekend worden
%
   if x < x1,
      yexact=x*(a1*x+b1)+c1+d1*exp(pe*(x-x1));
   else
      yexact=c2+d2*exp(pe*(x-1));
   end

