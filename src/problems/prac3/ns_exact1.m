%******************************************************************************
% DESCRIPTION:
%    This module computes the exact solution of 
%           dy/dx - k d2y/dx2 = max(2 - 5x,0)
%    on the interval 0 < x < 1, with boundary conditions y(0)=0, y(1)=1.
%    The solution is computed in the given grid point x.
%******************************************************************************
% INPUT:            x   coordinate grid point
%                   k   diffusioin coefficient
% OUTPUT:      yexact   exact solution in grid point
%******************************************************************************
%
function yexact=ns_exact(x,k)
%
pe=1/k;
% action to avoi underflow of exp(-pe)
if pe > 100,
   empe=0;
else
   empe=exp(-pe);
end
%
% general right-hand side is rl*(x1-x)
rl=5;
x1=0.4;
%
% on x < x1 is solution a1*x^2 + b1*x + c1 + d1*exp(pe*(x-x1))
% on x > x1 is solution c2 + d2*exp(pe*(x-1))
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
% solution can now be computed
%
   if x < x1,
      yexact=x*(a1*x+b1)+c1+d1*exp(pe*(x-x1));
   else
      yexact=c2+d2*exp(pe*(x-1));
   end

