function [x]=tdma(b,a,c,f,N);

% solve Ax=f using TDMA algorithm
%
% input: b,a,c = diagonals L,D,R of tridiagonal matrix
% input: f = r.h.s. vector
% input: N = dimension of problem
%
% output: x = solution vector

x=zeros(N,1);
alfa=zeros(N-1,1);
gamma=zeros(N-1,1);
y=zeros(N,1);

%step1: LU decomposition
alfa(1)=a(1);
gamma(1)=c(1)/alfa(1);
for i=2:N-1,
    alfa(i)=a(i)-b(i)*gamma(i-1);
    gamma(i)=c(i)/alfa(i);
end
alfa(N)=a(N)-b(N)*gamma(N-1);

%step2a: solve Ly=f
y(1)=f(1)/alfa(1);
for i=2:N,
    y(i)=(f(i)-b(i)*y(i-1))/alfa(i);
end

%step2b: solve Ux=y
x(N)=y(N);
for i=N-1:-1:1,
    x(i)=y(i)-gamma(i)*x(i+1);
end

