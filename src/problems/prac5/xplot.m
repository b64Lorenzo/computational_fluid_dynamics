function xplot(n,var,orient,plane)

% XPLOT 
%        XPLOT(N,VAR,ORIENT,PLANE) generates line plots of a number of 
%        variables along specified cross sections 
%        
%        N specifies the snapshot number
% 
%        VAR specifies which variable is plotted
%           1 - absolute velocity;
%           2 - horizontal velocity;
%           3 - vertical velocity;
%           4 - pressure;
%           5 - temperature (or label);
%           6 - absolute humidity;
%           7 - relative humidity;
%
%        ORIENT specifies whether the plane 
%           1 - horizontal
%           2 - vertical
%
%        PLANE specifies the index of the plane to be considered
%
%        Input to the routine are the files 'config.dat' and 'uvpf####.dat'

%        Filename: xplot.m
%        Date    : 5 November 1997
%        Author  : A.E.P. Veldman

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GLOBAL variables

global lq_snap lq_prt lq_xmode lq_xmin lq_xmax lq_xhold lq_legend nwetxt
global plotfg tel tekst t oudval

linevec=['-o','-*','-x','-+','-<','- '];
oud=get(lq_xhold,'Value');
if oud==0,
   clf reset
   tekst=[];
   tel=0;
%   plotmenu
else
   hold on
   tel=tel+1;
   if oudval==0
    if (round(t) <10)
     tekst1=['T= ', int2str(round(t))];
    elseif (round(t) <100)
     tekst1=['T=', int2str(round(t))];
    end
    tekst=cat(1,tekst,tekst1);
    legend(tekst,'NorthEastOutside')  
   end
end
oudval=oud;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read geometry parameters, labels and motion parameters

id=fopen('config.dat','r');
   geo=fscanf(id,'%e',[5,1]);
   icyl=geo(1); 
   xmin=geo(2);
   xmax=geo(3);
   ymin=geo(4);
   ymax=geo(5);
   nxy=fscanf(id,'%e',[2,1]);
   nx=nxy(1);
   ny=nxy(2);
% read cell centers
   xi=fscanf(id,'%e',[nx,1]);
   yj=fscanf(id,'%e',[ny,1]);
% read labels
   nf=fscanf(id,'%e',[nx,ny]);
% read motion parameters
   beweging=fscanf(id,'%e',[9,1]);
   ampl=beweging(1);
   freq=beweging(2);
   angle=beweging(3);
   omega=beweging(4);
   x0=beweging(5);
   y0=beweging(6);
if size(beweging,1) > 6,
   shks=beweging(7);
   stroke=beweging(8);
   L=beweging(9);
else
   shks=0;
end
fclose(id);

if xmin > 1.e-4,
% do not mirror picture
  icyl=0;
end  

% coordinates of cell edges

x = [2*xi(1)-xmin; 0.5*(xi(1:nx-1)+xi(2:nx)); 2*xi(nx)-xmax];
y = [2*yj(1)-ymin; 0.5*(yj(1:ny-1)+yj(2:ny)); 2*yj(ny)-ymax];

% mirror data when axisymmetric

imax=nx;
if icyl==1,
   xi = [-flipud(xi(2:nx)); xi(2:nx)];
   x = [-flipud(x(3:nx+1)); 0 ; x(3:nx+1)];
   xmin = - xmax;
   imax = 2*nx-2;
   nf = [ flipud(nf(2:nx,:)); nf(2:nx,:)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin plot loop

   number=sprintf('%04d', n);
   naam=['uvpf' number];
   humnaam=['humi' number];

% read u, v, p and f   

   id=fopen([naam '.dat'],'r');
      if id==-1,
         disp(['file ', naam,' not found']);
      end
      t=fscanf(id,'%e',[1,1]);
      uvpfdata=fscanf(id,'%e',[1,inf]);
      nuvpf=length(uvpfdata);
      if nuvpf==5*nx*ny,
         variables=5;
         uvpf=reshape(uvpfdata,5,nx*ny);
      elseif nuvpf==6*nx*ny,
         variables=6;
         uvpf=reshape(uvpfdata,6,nx*ny);
      elseif nuvpf==7*nx*ny,
         variables=7;
         uvpf=reshape(uvpfdata,7,nx*ny);
      elseif nuvpf==8*nx*ny,
         variables=8;
         uvpf=reshape(uvpfdata,8,nx*ny);
      else
         disp(['size of uvpf-file is incorrect']);
      end

      u=reshape(uvpf(1,:),nx,ny);
      v=reshape(uvpf(2,:),nx,ny);
      p=reshape(uvpf(3,:),nx,ny);
      f=reshape(uvpf(4,:),nx,ny);

      if variables==5,                    % temp./nf
         temp=reshape(uvpf(5,:),nx,ny);
      elseif variables==6,                % nf, temperature
         nfmov=reshape(uvpf(5,:),nx,ny);
         temp=reshape(uvpf(6,:),nx,ny);
      elseif variables==7,                % temp./nf, abs./rel. hum.
         temp=reshape(uvpf(5,:),nx,ny);
         abshum=reshape(uvpf(6,:),nx,ny);
         relhum=reshape(uvpf(7,:),nx,ny);
      elseif variables==8,                % nf, temp., abs./rel/ hum.
         nfmov=reshape(uvpf(5,:),nx,ny);
         temp=reshape(uvpf(6,:),nx,ny);
         abshum=reshape(uvpf(7,:),nx,ny);
         relhum=reshape(uvpf(8,:),nx,ny);
      end

      clear uvpf
   fclose(id);

% read absolute and relative humidity

   if var==6 | var==7,
      id=fopen([humnaam '.dat'],'r');
         if id==-1,
            disp(['file ', humnaam,' not found']);
         end
         hum=fscanf(id,'%e',[2,inf]);
         abshum=reshape(hum(1,:),nx,ny);
         relhum=100*reshape(hum(2,:),nx,ny);
         clear hum
      fclose(id);
   end

% mirror data when axisymmetric

   if icyl==1,
      u = [-flipud(u(2:nx,:)); u(2:nx,:)];
      v = [ flipud(v(2:nx,:)); v(2:nx,:)];
      p = [ flipud(p(2:nx,:)); p(2:nx,:)];
      f = [ flipud(f(2:nx,:)); f(2:nx,:)];
      temp = [ flipud(temp(2:nx,:)); temp(2:nx,:)];
      if var==6 | var==7,
	abshum = [ flipud(abshum(2:nx,:)); abshum(2:nx,:)];
	relhum = [ flipud(relhum(2:nx,:)); relhum(2:nx,:)];
      end;
   end

% show data only in liquid region

   q=sqrt(u.*u+v.*v);   
   vmax=max(max(q));
   leeg=find(f < 0.001);
   q(leeg)=NaN*ones(size(leeg));
   p(leeg)=NaN*ones(size(leeg));
   temp(leeg)=NaN*ones(size(leeg));
   if var==6 | var==7,
     abshum(leeg)=NaN*ones(size(leeg));
     relhum(leeg)=NaN*ones(size(leeg));
   end;

% make plot   

if orient==1
% horizontal cross section
 if (oud==1)
  telmod=mod(tel,6);
  if var==1
    if telmod==0
      plot(xi,sqrt(u(:,plane).^2+v(:,plane).^2),linevec(11:12))
    else
      plot(xi,sqrt(u(:,plane).^2+v(:,plane).^2),linevec(2*telmod-1:2*telmod))
    end
    txt='absolute velocity';
  elseif var==2
    if telmod==0
      plot(xi,u(:,plane),linevec(11:12))
    else
      plot(xi,u(:,plane),linevec(2*telmod-1:2*telmod))
    end
    txt='horizontal velocity';
  elseif var==3
    if telmod==0
      plot(xi,v(:,plane),linevec(11:12))
    else
      plot(xi,v(:,plane),linevec(2*telmod-1:2*telmod))
    end    
    txt='vertical velocity';
  elseif var==4
    if telmod==0
      plot(xi,p(:,plane),linevec(11:12))
    else
      plot(xi,p(:,plane),linevec(2*telmod-1:2*telmod))
    end    
    txt='pressure';    
  elseif var==5
    if mod(tel,size(a,2))==0
      plot(xi,temp(:,plane),linevec(11:12))
    else
      plot(xi,temp(:,plane),linevec(2*telmod-1:2*telmod))
    end    
    txt='temperature';    
  elseif var==6  
    if mod(tel,size(a,2))==0
      plot(xi,abshum(:,plane),linevec(11:12))
    else
      plot(xi,abshum(:,plane),linevec(2*telmod-1:2*telmod))
    end
    txt='absolute humidity';    
  elseif var==7
    if mod(tel,size(a,2))==0
      plot(xi,relhum(:,plane),linevec(11:12))
    else
      plot(xi,relhum(:,plane),linevec(2*telmod-1:2*telmod))
    end
    txt='relative humidity (%)';    
  end	
  if (get(lq_xmode,'Value') == 1)
    rmin=str2num(get(lq_xmin,'String'));
    rmax=str2num(get(lq_xmax,'String'));
    set(gca,'YLim',[rmin rmax]);
  else
    axis('auto')
  end
  xlabel('X','Fontsize',12);
  ylabel(txt,'Fontsize',12);
  title(['Y= ', num2str(yj(plane))],'Fontsize',16);
  if (round(t) <10)
     tekst1=['T= ', int2str(round(t))];
  elseif (round(t) <100)
     tekst1=['T=', int2str(round(t))];
  end
  tekst=cat(1,tekst,tekst1);
  legend(tekst,'NorthEastOutside')  
 else  
   if var==1
      plot(xi,sqrt(u(:,plane).^2+v(:,plane).^2))
      txt='absolute velocity';
   elseif var==2
      plot(xi,u(:,plane))
      txt='horizontal velocity';
   elseif var==3
      plot(xi,v(:,plane))
      txt='vertical velocity';
   elseif var==4
      plot(xi,p(:,plane))
      txt='pressure';
   elseif var==5
      plot(xi,temp(:,plane))
      txt='temperature';
   elseif var==6
      plot(xi,abshum(:,plane))
      txt='absolute humidity';
   elseif var==7
      plot(xi,relhum(:,plane))
      txt='relative humidity (%)';
   end
   if (get(lq_xmode,'Value') == 1)
      rmin=str2num(get(lq_xmin,'String'));
      rmax=str2num(get(lq_xmax,'String'));
      set(gca,'YLim',[rmin rmax]);
    else
      axis('auto')
    end
   xlabel('X','Fontsize',12);
   ylabel(txt,'Fontsize',12);
   title(['T= ',num2str(t), '     Y= ', num2str(yj(plane))],'Fontsize',16);
 end
else % orient==0, i.e. vertical cross section
 if (oud==1)
   telmod=mod(tel,6);
   if var==1
      if telmod==0
       plot(sqrt(u(plane,:).^2+v(plane,:).^2),yj,linevec(11:12))
      else
       plot(sqrt(u(plane,:).^2+v(plane,:).^2),yj,linevec(2*telmod-1:2*telmod))
      end     
      txt='absolute velocity';
   elseif var==2
      if telmod==0
        plot(u(plane,:),yj,linevec(11:12))
      else
        plot(u(plane,:),yj,linevec(2*telmod-1:2*telmod))
      end
      txt='horizontal velocity';
   elseif var==3
      if telmod==0
        plot(v(plane,:),yj,linevec(11:12))
      else
        plot(v(plane,:),yj,linevec(2*telmod-1:2*telmod))
      end      
      txt='vertical velocity';
   elseif var==4
      if telmod==0
        plot(p(plane,:),yj,linevec(11:12))
      else
        plot(p(plane,:),yj,linevec(2*telmod-1:2*telmod))
      end  
      txt='pressure';
   elseif var==5
      if mod(tel,size(a,2))==0
        plot(temp(plane,:),yj,linevec(11:12))
      else
        plot(temp(plane,:),yj,linevec(2*telmod-1:2*telmod))
      end    
      txt='temperature';
   elseif var==6
      plot(abshum(plane,:),yj)
      if mod(tel,size(a,2))==0
        plot(abshum(plane,:),yj,linevec(11:12))
      else
        plot(abshum(plane,:),yj,linevec(2*telmod-1:2*telmod))
      end
      txt='absolute humidity';
   elseif var==7
      if mod(tel,size(a,2))==0
        plot(relhum(plane,:),yj,linevec(11:12))
      else
        plot(relhum(plane,:),yj,linevec(2*telmod-1:2*telmod))
      end     
      txt='relative humidity';
   end
   if (get(lq_xmode,'Value') == 1)
      rmin=str2num(get(lq_xmin,'String'));
      rmax=str2num(get(lq_xmax,'String'));
      set(gca,'XLim',[rmin rmax]);
   else
      axis('auto')
   end
   xlabel(txt,'Fontsize',12);
   ylabel('Y','Fontsize',12);
   title(['X= ', num2str(xi(plane))],'Fontsize',16);
   if (round(t) <10)
     tekst1=['T= ', int2str(round(t))];
   elseif (round(t) <100)
     tekst1=['T=', int2str(round(t))];
   end
   tekst=cat(1,tekst,tekst1);
   legend(tekst,'NorthEastOutside')  
 else
% vertical cross section
   if var==1
      plot(sqrt(u(plane,:).^2+v(plane,:).^2),yj)
      txt='absolute velocity';
   elseif var==2
      plot(u(plane,:),yj)
      txt='horizontal velocity';
   elseif var==3
      plot(v(plane,:),yj)
      txt='vertical velocity';
   elseif var==4
      plot(p(plane,:),yj)
      txt='pressure';
   elseif var==5
      plot(temp(plane,:),yj)
      txt='temperature';
   elseif var==6
      plot(abshum(plane,:),yj)
      txt='absolute humidity';
   elseif var==7
      plot(relhum(plane,:),yj)
      txt='relative humidity';
   end
   if (get(lq_xmode,'Value') == 1)
      rmin=str2num(get(lq_xmin,'String'));
      rmax=str2num(get(lq_xmax,'String'));
      set(gca,'XLim',[rmin rmax]);
   else
      axis('auto')
   end
   xlabel(txt,'Fontsize',12);
   ylabel('Y','Fontsize',12);
   title(['T= ',num2str(t), '     X= ', num2str(xi(plane))],'Fontsize',16);
 end
end

% present figure

   drawnow
   

% define name of printfile
   
   if ~isempty(lq_prt),
      set(lq_prt,'String', [naam]);
   end

%zoom on

% end of plot loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
