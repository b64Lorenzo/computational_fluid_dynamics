figure(2)
subplot(211)
plot(-0.5,-0.5)
axis([1 tel 10^-6 10^1])
xlabel('iteration count')
ylabel('iteration error')
hold on
set(gca,'yscale','log')
semilogy(maxi)
if (maxi(tel-1)~=Inf)
text(1,5*10^-6,['   total number of iterations ', num2str(tel)]);
else
text(1,5*10^-6,'   the iterative method did not converge ');
end
hold off
%
%plot solutions
%
if rooster=='uni'
%
figure(2)
subplot(224)
set(0,'DefaultTextFontSize',10)
set(0,'DefaultAxesFontSize',10)
plot([0 x],[0 yrek],'r-',[0 x],[0;ydir],'k--',[0 x],[0 yex],'b-')
hold on
axis([0 1 0 1])
plot([0.03 0.05],[0.9 0.9],'b-')
text(0.1,0.9,'continuous solution')
plot([0.03 0.05],[0.8 0.8],'k-')
text(0.1,0.8,'direct solver')
plot([0.03 0.05],[0.7 0.7],'r-')
text(0.1,0.7,'iterative solver')
hold off
%
elseif rooster=='rek'
%
figure(2)
subplot(224)
set(0,'DefaultTextFontSize',10)
set(0,'DefaultAxesFontSize',10)
plot([0 xrek],[0 yrek],'r-',[0 xrek] ,[0;ydir],'k--',[0 xrek],[0 yexrek],'b-')
hold on
axis([0 1 0 1])
plot([0.03 0.05],[0.9 0.9],'b-')
text(0.1,0.9,'continuous solution')
plot([0.03 0.05],[0.8 0.8],'k--')
text(0.1,0.8,'direct solver')
plot([0.03 0.05],[0.7 0.7],'r-.')
text(0.1,0.7,'iterative solver')
hold off
%
end
