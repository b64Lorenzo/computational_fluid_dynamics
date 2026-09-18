%    PRACTICUM NUMERIEKE STROMINGSLEER
%
% OPGAVE 1    
% NAAM MODULE: cfd_2.m (oude naam: ns_1.m)
% DATUM:  september 2001 (oude versie: april 1994)
%**********************************************************************
% OMSCHRIJVING:
%    Deze module lost een convectie-diffusie vergelijking op
%    met diverse discretisatie methoden. De afzonderlijke methoden
%    zijn geimplementeerd in modules ns_1a.m, ns_1b.m, etc. 
%    De exacte oplossing wordt door ns_exact.m berekend
% INVOER:
% N     het aantal roosterpunten
% k     de diffusiecoefficient
% UITVOER:
% plaatjes van de oplossing voor acht methoden:
%  i) upwind; ii) smart upwind iii) B3; iv) QUICK; v) centraal; 
%  vi) upwind gerekt; vii) en viii) gerekt centraal
%***********************************************************************
clf
% oude arrays verwijderen
clear x yex yupw ylam xrek yrek yexrek ew A  
%
% invoer
%
N=input('give N:     ');
k=input('give k:     ');
%
Nex=1000;                      % Nex punten voor exacte oplossing
xex=1/Nex*[1:1:Nex];                          % definitie uniform rooster
for i=1:Nex
   yex(i)=ns_exact(xex(i),k);             % berekening exacte oplossing
end
%
x=1/N*[1:1:N];                  % definitie uniform rooster
yupw=zeros(1,N);                % initialiseren i.v.m. plots
ylam=zeros(1,N);
yrek=zeros(1,N);
%
% upwind
%
methode='U';
cfd_upw                      % berekening upwind oplossing - yupw
subplot(331)
plot(xex,yex,'b',x,yupw,'r')
axis([0 1 0 1])
title('upwind')
%
% smart-upwind
%
methode='S';
cfd_upw
subplot(332)
plot(xex,yex,'b',x,yupw,'r')
axis([0 1 0 1])
title('smart upwind')
%
% B3
%
lambda=0.5;
cfd_lam                     % berekening met lambda schema - ylam
subplot(334)
plot(xex,yex,'b',x,ylam,'r')
axis([0 1 0 1])
title('B3')
%
% quick
%
lambda=0.125;                 
cfd_lam                     
subplot(335)
plot(xex,yex,'b',x,ylam,'r')
axis([0 1 0 1])
title('QUICK')
%
% centraal
%
lambda=0.0;
cfd_lam
subplot(336)
plot(xex,yex,'b',x,ylam,'r')
axis([0 1 0 1])
title('central')
%
%pause
%clg
%
% gerekte rooster - exacte oplossing wordt in ns_1c.m bepaald
%
% upwind gerekt; Methode U
%
methode='U';
cfd_cen
subplot(337)
plot(xrek,yexrek,'b',xrek,yrek,'r')
axis([0 1 0 1])
title('upwind nonuniform')
%
% centraal gerekt; Methode A
%
methode='A';
cfd_cen
subplot(338)
plot(xrek,yexrek,'b',xrek,yrek,'r')
axis([0 1 0 1])
title('central nonuniform A')
%
% centraal gerekt; Methode B
%
methode='B';
cfd_cen
subplot(339)
plot(xrek,yexrek,'b',xrek,yrek,'r')
axis([0 1 0 1])
title('central nonuniform B')
%
[V,D]=eig(A);               % bepaling eigenwaarden coeff. matrix  
for i=1:N-1
   ew(i)=D(i,i);            % eigenwaarden opslaan in array ew 
end
%
%pause
%
% printen van de eigenwaarden
%
disp(['  ']);
disp(['eigenvalues of Method B']);
format long
ewB=ew'
format short

%legenda aanbrengen
subplot(333)
axis([0 1 0 1])
set(gca,'box','off')
set(gca,'ycolor',[1 1 1])
set(gca,'xcolor',[1 1 1])
set(gca,'ytick',[])
set(gca,'xtick',[])
text(0.2,0.85,['N=',num2str(N)])
text(0.2,0.65,['k=',num2str(k)])
hold on
plot([0.2,0.3],[0.35 0.35],'r-')
text(0.4,0.35,'discrete')
plot([0.2,0.3],[0.15 0.15],'b-')
text(0.4,0.15,'exact')
hold off

