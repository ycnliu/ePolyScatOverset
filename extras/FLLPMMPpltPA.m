%
% test of matlab
%
% Compute the polarization averaged distribution
%
clear;
fid=fopen('FLLPMMPPA.dat','r');
if fid==-1;
   disp('Could not find FLLPMMP.dat');
end;
np=40;
NPics=3;
thetandeg(1:NPics)=0; phindeg(1:NPics)=0;
thetandeg(1)=0.;  phindeg(1)=0.;
thetandeg(2)=90.; phindeg(2)=0.;
thetandeg(3)=90.; phindeg(3)=90.;
lmax=fscanf(fid,'%i', [1 1]);
 coefs = fscanf(fid, '%g', [6 inf]); 
 coefs=coefs';
disp('coefs read in ' );
disp(size(coefs));
theta= linspace(0,pi,np);
phi = linspace(0,2*pi, np);
[Phi,Theta] = meshgrid(phi,theta);
PHold =legendre(0,cos(theta))';
for i = 1:lmax;
   PHold = [PHold, legendre(i, cos(theta))'];
end;
dist=zeros(np,np);
for iPics =1:NPics;
thetan = thetandeg(iPics)*pi/180;
phin = phindeg(iPics)*pi/180;
PPHold = legendre(0, cos(thetan))';
for i = 1:2;
   PPHold = [PPHold, legendre(i, cos(thetan))'];
end;
for j=1:size(coefs,1);
   c = coefs(j, 1);
   l = coefs(j, 2);
   lp =coefs(j, 3);
   m = coefs(j, 4);
   mp = coefs(j, 5);
   sc = coefs(j, 6);
   mabs = abs(m);
   mpabs = abs(mp);
   ip = mabs+1+(l*(l+1))/2;
   ipp = mpabs+1+(lp*(lp+1))/2;
   if sc == 1;
      for k1 = 1:np;
         for k2 = 1:np;	
            dist(k1,k2)=dist(k1,k2)+c*PHold(k1, ip)*PPHold(1, ipp)*cos(m*phi(k2)+mp*phin);
         end;
      end;
   else;
      for k1 = 1:np;
         for k2 = 1:np;  
            dist(k1,k2)=dist(k1,k2)+c*PHold(k1, ip)*PPHold(1, ipp)*sin(m*phi(k2)+mp*phin);
         end;
      end;
   end;
end;
end;
dist=dist/3.;
mx = max(max(abs(dist)));
viewaz=[-40,0,0,90];
viewel=[20,90,0,0];
viewlab={'side','top','fy','fx'};
for iview = 1:length(viewaz);
titre=['PolAve',' ',viewlab{iview}];
titlefile='PolAve';
disp('Title');
disp(titre);
disp(titlefile);

% Use to only produce output file
% h=figure('visible','off');

% Use to produce output and also put on screen
h=figure('Name',titre,'NumberTitle','off')

construire_3D (Theta,Phi,dist,mx,titre,viewaz(iview),viewel(iview));
% titreblank=findstr(' ', titre);
% filestem=strcat(titre(1:titreblank(1)-1),titre(titreblank(1)+1:titreblank(2)-1));

      filestem='';
      for ii=1:length(titlefile);
         if titlefile(ii) ~= ' ';
           filestem = strcat(filestem, titlefile(ii));
         end;
         if length(filestem) >= 11;
            break;
         end;
      end;
      filestem = strcat(filestem,viewlab{iview})


% saveas(h,strcat(filestem,'.eps'),'epsc')
print(h,strcat(filestem,'.eps'),'-depsc','-r600');
% saveas(h,strcat(filestem,'.pdf'))
disp('Write files');
disp(strcat(filestem,'.eps'));
% disp(strcat(filestem,'.pdf'));
end;
% exit;


