%
%print eigenvalues 
%
figure(2)
subplot(223)
plot(-0.5,-0.5)    %commando to make sure the hold command works
axis([0 1 0 N])
set(gca,'box','off')
set(gca,'ycolor',[1 1 1])
set(gca,'xcolor',[1 1 1])
set(gca,'ytick',[])
set(gca,'xtick',[])
hold on
set(0,'DefaultTextFontSize',10)
text(0.1,N-0.5,'eigenvalues Jacobi matrix')
ewJacsort=sort(ewJac);
for i=1:N-1
% the eigenvalue 0
if   (abs(real(ewJacsort(i)))<10^(-10)) ...
   & (abs(imag(ewJacsort(i)))<10^(-10))
text(0.3,N-(i+1/2),['  ',num2str(0)])
% the pure imaginary eigenvalues
elseif (abs(real(ewJacsort(i)))<10^(-10)) ...
   & (imag(ewJacsort(i)) > 0)
text(0.3,N-(i+1/2),['  ',num2str(imag(ewJacsort(i))),'i'])
elseif (abs(real(ewJacsort(i)))<10^(-10)) ...
   & (imag(ewJacsort(i)) < 0)
text(0.3,N-(i+1/2),[num2str(imag(ewJacsort(i))),'i'])
% the real eigenvalues
elseif (abs(imag(ewJacsort(i)))<10^(-10)) ...
   & (real(ewJacsort(i)) > 0)
text(0.3,N-(i+1/2),['  ',num2str(real(ewJacsort(i)))])
elseif (abs(imag(ewJacsort(i)))<10^(-10)) ...
   & (real(ewJacsort(i)) < 0)
text(0.3,N-(i+1/2),num2str(real(ewJacsort(i))))
%
end
end
hold off
%