%Computation of the solution with the JOR method
%
fprintf('\n')
%
figure(2)
subplot(2,1,1)
axis([1 2 10^-6 10^1])
set(gca,'yscale','log')
xlabel('iteration count')
ylabel('iteration error')
hold on
%
tel=0;
k=0;
yrek(1:N-1)=0;
yrekold(1:N-1)=0;
maxi=1;
%
while (maxi > 1e-5)
   tel=tel+1;
   yrekold=yrek;
   yrek(1)=(-upp(1)*yrekold(2)+b(1))/dia(1);
   for i=2:N-2
     yrek(i)=(-upp(i)*yrekold(i+1)-low(i)*yrekold(i-1)+b(i))/dia(i); 
   end
   yrek(N-1)=(-low(N-1)*yrekold(N-2)+b(N-1))/dia(N-1);
   yrek=w*yrek+(1-w)*yrekold;
   maxi(tel)=max(max(abs(yrek-yrekold)));
% start writing data
   if tel>(k+1)*100 
    k0=100*k+1;
    k1=100*(k+1);
    axis([1 k1 10^-6 10^1])
    set(gca,'yscale','log')
    semilogy([k0:k1],maxi([k0:k1]))
    pause(0.01)
    k=k+1;
   end
   if rem(tel,5)==0
    fprintf('n %4.0f      error %3.8e \n',tel,maxi(tel))
   end
% end writing data
end
%
hold off
fprintf('total number of iterations  %4.0f  ',tel)
fprintf('\n')
fprintf('\n')
fprintf('\n')
















