function [x]=pdma(b2,b1,a,c1,c2,f,N);

% solve Ax=f using PDMA algorithm
%
% input: b2,b1,a,c1,c2 = diagonals of penta-diagonal matrix
%        a=diag(A);
%        c1=[diag(A,1);0];
%        c2=[diag(A,2);0;0];
%        b1=[0;diag(A,-1)];
%        b2=[0;0;diag(A,-2)];
% input: f = r.h.s. vector
% input: N = dimension of problem
%
% output: x = solution vector

x=zeros(N,1);

alfa=zeros(N,1);
gamma=zeros(N-1,1);
delta=zeros(N-2,1);
beta=zeros(N,1);
y=zeros(N,1);

%step1: LU decomposition
alfa(1)=a(1);
gamma(1)=c1(1)/alfa(1);
delta(1)=c2(1)/alfa(1);
beta(2)=b1(2);
alfa(2)=a(2)-beta(2)*gamma(1);
gamma(2)=( c1(2)-beta(2)*delta(1) )/alfa(2);
delta(2)=c2(2)/alfa(2);
for k=3:N-2
    beta(k)=b1(k)-b2(k)*gamma(k-2);
    alfa(k)=a(k)-b2(k)*delta(k-2)-beta(k)*gamma(k-1);
    gamma(k)=( c1(k)-beta(k)*delta(k-1) )/alfa(k);
    delta(k)=c2(k)/alfa(k);
end

beta(N-1)=b1(N-1)-b2(N-1)*gamma(N-3);
alfa(N-1)=a(N-1)-b2(N-1)*delta(N-3)-beta(N-1)*gamma(N-2);
gamma(N-1)=( c1(N-1)-beta(N-1)*delta(N-2) )/alfa(N-1);
beta(N)=b1(N)-b2(N)*gamma(N-2);
alfa(N)=a(N)-b2(N)*delta(N-2)-beta(N)*gamma(N-1);

%step2a: solve Ly=f
y(1)=f(1)/alfa(1);
y(2)=(f(2)-beta(2)*y(1))/alfa(2);
for k=3:N
    y(k)=( f(k)-b2(k)*y(k-2)-beta(k)*y(k-1) )/alfa(k);
end

%step2b: solve Ux=y
x(N)=y(N);
x(N-1)=y(N-1)-gamma(N-1)*x(N);
for k=N-2:-1:1
    x(k)=y(k)-gamma(k)*x(k+1)-delta(k)*x(k+2);
end

