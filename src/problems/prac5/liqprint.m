function liqprint(or,fmt)

% LIQPRINT prints current figure 
%         LIQPRINT(OR,FMT) exports the current figure to an epsc-file.
%         The printfile is named by the process that created the figure, 
%         but it can be altered manually in the Liquid Control Menu (the 
%         corresponding control handle is LQ_PRT).
%
%         OR controls the way in which the geometry is plotted:
%              OR = 1: portrait,  2: landscape
%         FMT control the format of the picture-file:
%              FMT = 1: epsc,  2: jpg (gif with Matlab version 4)

%         Filename: liqprint.m
%         Date    : 27 March 1999
%         Author  : A.E.P. Veldman

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GLOBAL variables

global lq_prt 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
naam=get(lq_prt,'String');


if or == 1,
   orient portrait
else
   orient landscape
end

if fmt == 1,
   format = ' -depsc ';
   ext = '.epsc';
elseif fmt == 2,
   vernum=version;
   if vernum(1)=='4',
      format = ' -dgif8 ';
      ext = '.gif';
   else
      format = ' -djpeg ';
      ext = '.jpg';
   end
else
   disp(['wrong file format']);
end

disp(['the print-file  ', naam, ext, '  is being made -- please wait'])
eval(['print ', format, naam, ext])
disp(['print ', format, naam, ext '  has been finished'])
