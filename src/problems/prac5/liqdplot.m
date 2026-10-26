
function liqdplot(nb,ne,icol,imovie)

% LIQDPLOT 
%        LIQDPLOT(NB,NE,ICOL,IMOVIE) plots the liquid region in various 
%        ways at specified instants in time, defined by the interval [NB,NE].
%        
%        ICOL specifies how the liquid region is colored: 
%           1 - vof-function; 
%           2 - absolute velocity;
%           3 - horizontal velocity;
%           4 - vertical velocity;
%           5 - pressure;
%           6 - vorticity;
%           7 - temperature (or label);
%           8 - absolute humidity;
%           9 - relative humidity;
%          10 - fog indicator.
%
%        The color- and contour range can be manually prescribed,
%        as well as the scaling of velocity vectors.
%
%        The individual plots can be combined into a Matlab movie or be 
%        transformed into an Mpeg file. Making of an Mpeg file requires 
%        the availability of the Matlab routine mpgwrite.mexhp7 (or its 
%        equivalent on other platforms). 
%
%        IMOVIE specifies whether a screen movie should be made (IMOVIE > 0),
%           -1 - live mode;
%            2 - save as a Matlab-file;
%            3 - export as an mpeg movie;
%            4 - export as an AVI movie;
%            5 - export as jpg files.
%
%        If available streamfunction and vorticity data are presented.
%
%        All plots can be exported as an *.epsc file.
%
%        Moving boundaries automatic when input file gives label 
%
%        Input to the routine are the files 'config.dat', 'uvpf####.dat'
%        and 'psiomega.dat'.
%
%        Automatic adjustment to Matlab versions; in particular the major graphics step from R2014a to R2014b.

%        Filename: liqdplot.m
%        Date    : 5 March 2001
%        Last change: January 2016
%        Author  : A.E.P. Veldman

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axiszoom=[];
if (~isempty(get(gcf,'CurrentAxes'))) axiszoom=axis; end;
clf reset

% GLOBAl variables

global lq_snap lq_prt lq_mpg lq_quiver lq_relscale lq_absscale
global lq_cmode lq_cmin lq_cmax lq_live lq_sel lq_zoom
global lq_stride lq_bar lq_cntr lq_clbl x y data M
global lq_1axismin lq_1axismax lq_2axismin lq_2axismax lq_keepaspect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read geometry parameters, labels and motion parameters

% If in live mode, wait until config.dat exists

if (imovie==-1),
  id=fopen('config.dat');
  counter=1;
  while (id==-1 & counter<60),
    pause(1);
    id=fopen('config.dat');
    counter=counter+1;
    if strcmp(get(lq_live, 'String'), 'LIVE'), counter=60; end;
  end;
  if (id~=-1),
    fclose(id);
  else
    return;
  end;
end;

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

% define axis

if shks==0,
   scale=1.0;
   xamp=scale*(x0+ampl*cos(angle));
   yamp=scale*(y0+ampl*sin(angle));
else
   scale=0.2;
   xamp=0.0;
   yamp=scale*0.5*stroke;
end

set(gca,'Position', [0.04 0.04 0.92 0.87]);
vernum=version;
% default plot area
xlim=[x(1)-xamp x(imax+1)+xamp];
ylim=[y(1)-yamp y(ny+1)+yamp];
% custom plot area
eval('zoomval=get(lq_zoom, ''Value'');', 'zoomval=1;');
if isempty(zoomval), zoomval=1; end;
if (zoomval==2),
  xlim=[str2num(get(lq_1axismin, 'String')), ...
	str2num(get(lq_1axismax, 'String'))];
elseif (zoomval==3),
  ylim=[str2num(get(lq_1axismin, 'String')), ...
	str2num(get(lq_1axismax, 'String'))];
elseif (zoomval==4),
  xlim=[str2num(get(lq_1axismin, 'String')), ...
	str2num(get(lq_1axismax, 'String'))];
  ylim=[str2num(get(lq_2axismin, 'String')), ...
	str2num(get(lq_2axismax, 'String'))];
elseif (zoomval==5),
  if (~isempty(axiszoom)), xlim=axiszoom(1:2); ylim=axiszoom(3:4); end;
end;

eval('aspect=get(lq_keepaspect, ''Value'');', 'aspect=1;');
if isempty(aspect), aspect=1; end;
if vernum(1)=='4',
   axis([xlim ylim]);
   if (aspect), axis('equal'); end;
   axis('off');
else  
  set(gca,'Visible','off', 'XLim', xlim, 'YLim', ylim);
  if (aspect), 
    set(gca,'DataAspectRatio',[1 1 1]); 
  else
    set(gca,'DataAspectRatioMode', 'auto'); 
  end;
end;

% define colormap

colormap(copper);
map=colormap;
bluemap=[map(:,3) map(:,2) map(:,1)];

% initialize memory for movie

if imovie >= 1,
   M=moviein(ne-nb+1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot loop

for n=nb:ne,

% If in live mode, check whether the next plot (i.e. n+1) is
% already present. In this way it is not possible that a file is
% being caught while it is still being created.
   if (imovie==-1),
     name=['uvpf',sprintf('%04d', n+1),'.dat'];
     id=fopen(name);
     counter=1;
     while (id==-1 & counter < 60),
       pause(0.5);
       id=fopen(name);
       counter=counter+1;
       if strcmp(get(lq_live, 'String'), 'LIVE'), counter=60; end;
     end;
     if (id==-1), % file does not exist
       disp(['file ', name, ' not found']);
       return;
     else         % file exists, continue with plot
       fclose(id);
       set(lq_sel, 'String', num2str(n));
     end;
     if strcmp(get(lq_live, 'String'), 'LIVE'), return; end;
   end;
       
   number=sprintf('%04d', n);
   naam=['uvpf' number];
   humnaam=['humi' number];

% read u, v, p and f and other variables

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

% compute vorticity  

   vort=NaN*ones(size(nf));
   emp = 0.001;
   if icol==6,
      for j=2:ny-1,
      for i=2:nx-1,
         if ( f(i+1,j) > emp & f(i-1,j) > emp & ...
                   f(i,j+1) > emp & f(i,j-1) > emp ),
            vort(i,j)=-(u(i,j+1)-u(i,j-1))/(yj(j+1)-yj(j-1)) ...
                   + (v(i+1,j)-v(i-1,j))/(xi(i+1)-xi(i-1));
         end
      end
      end
   end

% read absolute and relative humidity

   if (icol>=8 & icol<=10),
     if (variables~=7 & variables~=8),
       id=fopen([humnaam '.dat'],'r');
       if id==-1,
	 disp(['file ', humnaam,' not found']);
       end
       hum=fscanf(id,'%e',[2,inf]);
       fclose(id);
       abshum=reshape(hum(1,:),nx,ny);
       relhum=reshape(hum(2,:),nx,ny);
     end;
     fog=(relhum-1).*(relhum>1.0);
     relhum=100*relhum;
     clear hum
   end;

% mirror data when axisymmetric

   if icyl==1,
      u = [-flipud(u(2:nx,:)); u(2:nx,:)];
      v = [ flipud(v(2:nx,:)); v(2:nx,:)];
      p = [ flipud(p(2:nx,:)); p(2:nx,:)];
      f = [ flipud(f(2:nx,:)); f(2:nx,:)];
      temp = [ flipud(temp(2:nx,:)); temp(2:nx,:)];
      vort = [ -flipud(vort(2:nx,:)); vort(2:nx,:)];
      if (variables==6 | variables==8),
         nfmov = [ flipud(nfmov(2:nx,:)); nfmov(2:nx,:)];
      end
      if (icol>=8 & icol<=10),
         abshum = [ flipud(abshum(2:nx,:)); abshum(2:nx,:)];
         relhum = [ flipud(relhum(2:nx,:)); relhum(2:nx,:)];
         fog    = [ flipud(fog(2:nx,:)); fog(2:nx,:)];
      end
   end

   
% treatment of moving boundaries

   if (variables==5 | variables==7),
%     check whether input describes `nf' or `temp' (nf is always integer)
      rnf=round(temp);
      if rnf == temp,
         nf=rnf;
      end
   elseif (variables==6 | variables==8),
% new `nf' is available
      nf=nfmov;
   end


% show data only in liquid region

   q=sqrt(u.*u+v.*v);
   vmax=max(max(q));
   leeg=find(f < 0.001);
   q(leeg)=NaN*ones(size(leeg));
   p(leeg)=NaN*ones(size(leeg));
   temp(leeg)=NaN*ones(size(leeg));
   if (icol>=8 & icol<=10),
      abshum(leeg)=NaN*ones(size(leeg));
      relhum(leeg)=NaN*ones(size(leeg));
      fog(leeg)=NaN*ones(size(leeg));
   end

% convert motion to position on screen
   
   xoffset=scale*(x0+ampl*cos(angle))*cos(freq*t);
   xiplot=xi+xoffset;
   xplot=x+xoffset;

  if shks==0,
    yoffset=scale*(y0*sin(freq*t)+ampl*sin(angle)*cos(freq*t));
  else
    r=0.5*stroke;
    lambda=r/L;
    w=shks*pi/30;
    fi=w*t+pi;
    yoffset=-r*cos(fi) + L*(sqrt(1.0-lambda^2*(sin(fi))^2) - 1.0);
    yoffset=scale*yoffset;
  end
   yjplot=yj+yoffset;
   yplot=y+yoffset;

% make plot   

   caxis('auto')
   cmanual = get(lq_cmode,'Value');
   if cmanual == 1,                        
      cmin=str2num(get(lq_cmin,'String'));
      cmax=str2num(get(lq_cmax,'String'));
   end

   icntr = get(lq_cntr,'Value') - 1;

   data=[];
   if icntr==0,                                 % colorplot
      if icol==1,
         pcolor(xiplot,yjplot,f')
         colormap(bluemap)
	 data=f;
      elseif icol==2,
         pcolor(xiplot,yjplot,q')
         colormap(jet)
	 data=q;
      elseif icol==3,
         pcolor(xiplot,yjplot,u')
         colormap(jet)
	 data=u;
      elseif icol==4,
         pcolor(xiplot,yjplot,v')
         colormap(jet)
	 data=v;
      elseif icol==5,
         pcolor(xiplot,yjplot,p')
         colormap(jet)
	 data=p;
      elseif icol==6,
         pcolor(xiplot,yjplot,vort')
         colormap(jet)
	 data=vort;
      elseif icol==7,
         pcolor(xiplot,yjplot,temp')
         colormap(jet)
	 data=temp;
      elseif icol==8,
	 pcolor(xiplot,yjplot,abshum')
	 colormap(jet)
	 data=abshum;
      elseif icol==9 ,
	 pcolor(xiplot,yjplot,relhum')
	 colormap(jet)
	 data=relhum;
      elseif icol==10,
	 pcolor(xiplot,yjplot,fog')
	 colormap(1-gray)
	 caxis([0 1])
	 data=fog;
      end
      if cmanual == 1,                         % manual coloring
         caxis([cmin cmax]);
      end
   elseif (icntr==1 & cmanual==0),             % contourplot auto
      if icol==1,
         cc=contour(xiplot,yjplot,f');
	 data=f;
      elseif icol==2,
         cc=contour(xiplot,yjplot,q');
	 data=q;
      elseif icol==3,
         cc=contour(xiplot,yjplot,u');
	 data=u;
      elseif icol==4,
         cc=contour(xiplot,yjplot,v');
	 data=v;
      elseif icol==5,
         cc=contour(xiplot,yjplot,p');
	 data=p;
      elseif icol==6,
         cc=contour(xiplot,yjplot,vort');
	 data=vort;
      elseif icol==7,
         cc=contour(xiplot,yjplot,temp');
	 data=temp;
      elseif icol==8,
         cc=contour(xiplot,yjplot,abshum');
	 data=abshum;
      elseif icol==9 ,
         cc=contour(xiplot,yjplot,relhum');
	 data=relhum;
      elseif icol==10,
         cc=contour(xiplot,yjplot,fog');
	 data=fog;
      end
   elseif (icntr==1 & cmanual==1),            % contourplot manual 
      V= cmin + [0:10]*(cmax-cmin)/10;
      if icol==1,
         cc=contour(xiplot,yjplot,f',V);
	 data=f;
      elseif icol==2,
         cc=contour(xiplot,yjplot,q',V);
	 data=q;
      elseif icol==3,
         cc=contour(xiplot,yjplot,u',V);
	 data=u;
      elseif icol==4,
         cc=contour(xiplot,yjplot,v',V);
	 data=v;
      elseif icol==5,
         cc=contour(xiplot,yjplot,p',V);
	 data=p;
      elseif icol==6,
         cc=contour(xiplot,yjplot,vort',V);
	 data=vort;
      elseif icol==7, 
         cc=contour(xiplot,yjplot,temp',V);
	 data=temp;
      elseif icol==8,
         cc=contour(xiplot,yjplot,abshum',V);
	 data=abshum;
      elseif icol==9 ,
         cc=contour(xiplot,yjplot,relhum',V);
	 data=relhum;
      elseif icol==10,
         cc=contour(xiplot,yjplot,fog',V);
	 data=fog;
      end
   end

   shading interp

   if icol==1,
      titel='VOF';
   elseif icol==2,
      titel='velocity';
   elseif icol==3,
      titel='horizontal velocity';
   elseif icol==4,
      titel='vertical velocity';
   elseif icol==5,
      titel='pressure';
   elseif icol==6,
      titel='vorticity';
   elseif icol==7,
      titel='streamfunction';
   elseif icol==8,
      titel='absulute humidity';
   elseif icol==9 ,
      titel='relative humidity (%)';
   elseif icol==10,
      titel='fog indicator';
   end
   
   if vernum(1)=='4',  % version Matlab 4 
      axis([xlim ylim])
      if (aspect), axis('equal'); end;
      axis('off')
      title(['t = ',sprintf('%6.3e',t),' : ',titel],'Fontsize',16);
   else               % later Matlab version > 4
      set(gca, 'Visible','off', 'XLim', xlim, 'YLim', ylim, ...
          'Title',text('String',['t = ',sprintf('%6.3e',t),' : ',titel], ...
          'Fontsize',16))
% Jan 2016: changes in graphics handles for R2014b (= 8.4.0) and later
      if ~verLessThan('matlab','8.4.0')   
         ax=gca;
         ax.Title.Visible='on';
         ax.Title.FontWeight='normal';
         ax.Title.FontSize=12; 
      end
% end R2014b
      if (aspect), 
	set(gca,'DataAspectRatio',[1 1 1]); 
      else
	set(gca,'DataAspectRatioMode', 'auto'); 
      end;
   end

   
   hoek=0;
%   hoek=10*n;
   view([hoek 90])
   
   cb=[];
   ibar=get(lq_bar,'Value');               % show colorbar
   if imovie==0 & icol>1 & icntr==0 & ibar==1,
      if ~(icol==10),
	cb=colorbar;
      end
   end

   ilabel=get(lq_clbl,'Value');            % label contour values
   if imovie==0 & icntr==1 & ilabel==1,
      clabel(cc)
   end

   hold on

% draw geometry

   for j=1:ny,
      for i=1:imax,
         xl=xplot(i);
         xr=xplot(i+1);
         yb=yplot(j);
         yt=yplot(j+1);
         if (nf(i,j)==9) 
            fill([xl xr xr xl],[yb yb yt yt],[0.55 0.35 0.25] )
         elseif (nf(i,j)==8)
            fill([xl xr xr xl],[yb yb yt yt],'y')
         elseif (nf(i,j)==7)
%            fill([xl xr xr xl],[yb yb yt yt],[0.0 0.0 0.0] )
         end
      end
   end       

% add velocity vectors

   if get(lq_quiver, 'Value')>1,
      scale=1.5;
      stride=1;
      if (get(lq_quiver,'Value')==3)
         scale=get(lq_relscale,'Value');
         scale=1.5*exp(scale);
         stride=get(lq_stride,'Value');
      end
      if (get(lq_quiver,'Value')==4)
         scale=get(lq_absscale,'String');
	 scale=str2num(scale);
         stride=get(lq_stride,'Value');
      end
      [mx,my]=size(u);
      for j=1:floor(my/stride)
         yjstride(j)=yjplot(1+stride*(j-1));
         for i=1:floor(mx/stride)
            ustride(i,j)=u(1+stride*(i-1),1+stride*(j-1));
            vstride(i,j)=v(1+stride*(i-1),1+stride*(j-1));
            xistride(i)=xiplot(1+stride*(i-1));
         end
      end      
      hold on
      if (get(lq_quiver,'Value')==4)
	quiver(xistride,yjstride,scale*ustride',scale*vstride',0,'b')
      else
	quiver(xistride,yjstride,ustride',vstride',scale*vmax,'b')
      end
   end

% present figure

%  if (get(lq_keepzoom, 'Value')==1 & ~isempty(axiszoom)),
%    axis(axiszoom);
%  end;
  drawnow

% there seems to be a bug in Matlab6; if the values of the colorbar
% are scaled (by a 10-power), then the x10^? is not displayed. This
% statement should help.
  if (~isempty(cb) & vernum(1)=='6'), 
    set(cb, 'YTickLabelMode', 'auto');
  end;
   
% make frame of movie

   if imovie >=1,
     if vernum(1)=='6',
       M(:,n-nb+1)=getframe(gcf); %the gcf parameter is necessary for matlab6
     else  
       M(:,n-nb+1)=getframe; 
     end
   end

% release figure   

   hold off

% define name of printfile
   
   if ~isempty(lq_prt),
      set(lq_prt,'String', [naam]);
   end

% save files of individual frames   
   if imovie==5,
      if ~isempty(lq_mpg),
         filmnaam=get(lq_mpg,'String');
      else
         filmnaam='uvpf';
      end

      if vernum(1)=='4',
         disp(['frame is saved as ', filmnaam,number, '.gif']);
         eval(['print -dgif8 ', filmnaam, number,'.gif']);
         disp(['frame ', number, ' ready']);
      else
         disp(['frame is saved as ', filmnaam,number, '.jpg']);
         eval(['print -djpeg ', filmnaam, number,'.jpg']);
         disp(['frame ', number, ' ready']);
      end
   end

% save current plot number to file 'veilig' in live mode
   if imovie==-1,
     eval(['!echo ', num2str(n), ' >veilig']);
   end;
   
end % end of plot loop


if imovie==0 & hoek==0,
   zoom on
end
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% show and export movie

if imovie >=1,
  if vernum(1)=='6',
    movie(gcf,M,3,5); %the gcf parameter is necessary for matlab6
  else
    movie(M,3,5);
  end
end

if imovie==2,

% save movie as matlab-file

   if ~isempty(lq_mpg),
     matfile=get(lq_mpg,'String');
   else
     matfile='matfile';
   end

   eval(['save ', matfile, ' M']);
   disp(['movie saved as file  ', matfile, '.mat']); 

elseif imovie==3,

% save movie as mpg-file

     disp('an mpeg movie is now being made -- please wait');

     if ~isempty(lq_mpg),
       mpgfile=get(lq_mpg,'String');
     else
       mpgfile='mpgfile';
     end

     if icol==1,
       eval(['mpgwrite(M,bluemap,mpgfile);']);
     elseif icol<=9
       eval(['mpgwrite(M,jet,mpgfile);']);
     else
       eval(['mpgwrite(M,1-gray,mpgfile);']);
     end

     disp(['mpeg movie  ', mpgfile,'.mpg  ready'])

elseif imovie==4,

% save movie as avi-file

     disp('an AVI movie is now being made -- please wait');

     if ~isempty(lq_mpg),
       mpgfile=get(lq_mpg,'String');
     else
       mpgfile='avifile';
     end

     if icol==1,
       eval(['movie2avi(M,mpgfile,''colormap'',bluemap);']);
     elseif icol<=9
       eval(['movie2avi(M,mpgfile,''colormap'',jet);']);
     else
       eval(['movie2avi(M,mpgfile,''colormap'',1-gray);']);
     end

     disp(['AVI movie  ', mpgfile,'.avi  ready'])

end  % end of movie-making

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



