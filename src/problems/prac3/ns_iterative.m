delete(get(figure(2),'children'))
%
% input
%
N=10;
k=str2num(kstr);
w=str2num(wstr);
%
if (rooster=='uni')
%
x=1/N*[1:1:N];                          % definition uniform grid
for i=1:N
   yex(i)=ns_exact(x(i),k);             % computation exact solution
end
%
h=1/N;                         % mesh size uniform grid
%
if methode=='Uuni'
%
% computation of coefficients, upwind
%
   upp=-k/(h^2)*ones(1,N-1);               % upper diagonal
   dia=(2*k/(h^2)+1/h)*ones(1,N-1);        % diagonal
   low=(-k/(h^2)-1/h)*ones(1,N-1);         % lower diagonal
%
elseif methode=='Cuni'
%
% computation of coefficients, central
%
   lambda=0.0;
   upp=(0.5/h-lambda/h-k/(h^2))*ones(1,N-1);
   dia=(2*k/(h^2)+3*lambda/h)*ones(1,N-1);
   low=(-0.5/h-3*lambda/h-k/(h^2))*ones(1,N-1);
%
% first equation central
%
   upp(1)=0.5/h-k/(h^2);
   dia(1)=2*k/(h^2);
%
end
%
% right-hand side
%
M=round(0.4/h);
b=2-5*h*(1:1:M);
b(M+1:N-1)=zeros(1,N-M-1);
b(N-1)=b(N-1)-upp(N-1);         % boundary condition on right
%
elseif (rooster=='rek')
%
% computation of stretched grid; N/2 meshes left and right of x = 1-Lk
% good choice L=5
%
L=5;
%
hcrs=2.*(1.-L*k)/N;           % coarse meshes left
hfin=2*L*k/N;                 % fine meshes right
for i=1:N/2
   xrek(i)=i*hcrs;            % grid points left
end
for i=N/2+1:N 
   xrek(i)=1.-(N-i)*hfin;     % grid points right
end
% 
% exact solution on stretched grid
%
for i=1:N
yexrek(i)=ns_exact(xrek(i),k);
end
%
if methode=='Urek'
%
% computation coefficients, upwind
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
%
else
%
% computation coefficients, central
%
   upp(1:N/2-1)=(0.5/hcrs-k/(hcrs^2))*ones(1,N/2-1);
   dia(1:N/2-1)=2.*k/(hcrs^2)*ones(1,N/2-1);
   low(1:N/2-1)=-0.5/hcrs-k/(hcrs^2)*ones(1,N/2-1);
%
   i=N/2;
   if methode=='Arek' 
      upp(i)=1./(hcrs+hfin)-2.*k/(hfin*(hfin+hcrs));
      dia(i)=2.*k/(hfin*hcrs);
      low(i)=-1./(hcrs+hfin)-2.*k/(hcrs*(hcrs+hfin));
   end
   if methode=='Brek' 
      upp(i)=(hcrs-2.*k)/(hfin*(hfin+hcrs));
      dia(i)=(hfin-hcrs+2.*k)/(hfin*hcrs);
      low(i)=(-hfin-2.*k)/(hcrs*(hcrs+hfin));
   end
%
   upp(N/2+1:N-1)=(0.5/hfin-k/(hfin^2))*ones(1,N/2-1);
   dia(N/2+1:N-1)=2.*k/(hfin^2)*ones(1,N/2-1);
   low(N/2+1:N-1)=(-0.5/hfin-k/(hfin^2))*ones(1,N/2-1);
%
end
%
% right-hand side
%
M=round(0.4/hcrs);
b=2.-5.*hcrs*(1:1:M);
b(M+1:N-1)=zeros(1,N-M-1);
b(N-1)=b(N-1)-upp(N-1);       % boundary condition right
%
end
%
%computation eigenvalues of Jacobi matrix
%
D=diag(dia(1:N-1));
LpU=diag(-upp(1:N-2),1)+diag(-low(2:N-1),-1);
A=diag(upp(1:N-2),1)+diag(dia(1:N-1),0)+diag(low(2:N-1),-1);
Jac=inv(D)*LpU;
[evmat,ewmat]=eig(Jac);
for i=1:N-1
   ewJac(i)=ewmat(i,i);   % save eigenvalues in array ewJac
end
%
plotev
%
% compute solution with iterative method
%
if itmethode=='JOR'
ns_jor
elseif itmethode=='SOR'
ns_sor
end
yrek(N)=1;
%
% compute solution with direct method
%
ydir=A\b';
ydir(N)=1;
%
plotsol
%


