function yy = HLLPMMPpltva(inpfile)
%
% test of matlab
%
fid=fopen(inpfile,'r');
if fid==-1;
   disp('Could not find file '); disp(inpfile);
end;
np=40
NPics=fscanf(fid,'%i', 1);
nvax=fscanf(fid,'%i', 1);
thetandeg(1:NPics)=0; phindeg(1:NPics)=0;
mu0(1:NPics)=0; ipSel(1:NPics)=0; itSel(1:NPics)=0;
fixphase(1:NPics)=0;
for i = 1:NPics
   angsn=fscanf(fid,'%g', [6, 1]);
   thetandeg(i) = angsn(1);   % polar direction of the field
   phindeg(i) = angsn(2);     % azimuthal direction of the field
   mu0(i) = angsn(3);         % type of field mu0 = +/- 1 for circularly polarized
                              %                   = 0 linearly polarized
   
   ipSel(i) = angsn(4);       % =1 for length, =2 for velocity
   itSel(i) = angsn(5);       % = the target component to use of a degenerate set
                              % if there is no degeneracy then this should be =1
   fixphase(i) = angsn(6);    % phase angle in degrees to adjust these matrix elements by exp(i*phase)
  blank = fgetl(fid);
  blank = fgetl(fid);
  if i==1;
     label=blank;
  else;
     label = strvcat(label, blank);
  end;
end;
lmax=fscanf(fid,'%i', [1 1]);
 coefs = fscanf(fid, '%g', [7 inf]); 
 coefs=coefs';
disp('coefs read in ' );
disp(size(coefs));
disp(coefs(1,:));
theta= linspace(0,pi,np);
phi = linspace(0,2*pi, np);
[Phi,Theta] = meshgrid(phi,theta);
PHold =legendre(0,cos(theta))';
for i = 1:lmax;
   PHold = [PHold, legendre(i, cos(theta))'];
end;
ij=0;
for l=0:lmax;
   lfactor=sqrt((2*l+1)/(4.*pi));
   for m=0:l;
      ij=ij+1;
      if m > 0;
         lfactor = lfactor*sqrt(1/((l+m)*(l+1-m)));
      end;
      PHold(:,ij) = lfactor*PHold(:,ij);
   end;
end;

for iPics =1:NPics
   thetan = thetandeg(iPics)*pi/180;
   phin = phindeg(iPics)*pi/180;
   PPHold = complex([0 0 0]);
   if mu0(iPics) == -1;   % mu0=-1 and mu0=1 are reversed since our matrix elements
                          % are backwards
      PPHold(3) = 0.5*(1+cos(thetan))*complex(cos(phin), -sin(phin)); % m=1
      PPHold(2) = 0.5*sqrt(2.0)*sin(thetan);                          % m=0
      PPHold(1) = 0.5*(1-cos(thetan))*complex(cos(phin), sin(phin));  % m=-1
   elseif mu0(iPics) == 0;
      PPHold(3) = -0.5*sqrt(2.0)*sin(thetan)*complex(cos(phin), -sin(phin));   % m=1
      PPHold(2) = cos(thetan);                                        % m=0
      PPHold(1) = 0.5*sqrt(2.0)*sin(thetan)*complex(cos(phin), sin(phin));  % m=-1
   else;    % mu0(iPics) = 1
      PPHold(3) = 0.5*(1-cos(thetan))*complex(cos(phin), -sin(phin)); % m=1
      PPHold(2) = -0.5*sqrt(2.0)*sin(thetan);                         % m=0
      PPHold(1) = 0.5*(1+cos(thetan))*complex(cos(phin), sin(phin));  % m=-1
   end;

   dist=zeros(np,np);
   phase=zeros(np,np);
   amplitude=complex(zeros(np,np));
   for j=1:size(coefs,1);
      m = coefs(j, 1);
      l = coefs(j, 2);
      mu = coefs(j, 3);
      ip = coefs(j, 4);
      it = coefs(j, 5);
      c = complex(coefs(j, 6), coefs(j, 7));
      mabs = abs(m);

      if ipSel(iPics) == ip & itSel(iPics) == it;
         ip = mabs+1+(l*(l+1))/2;
         ipp = mu+2;
         for k1 = 1:np;
            for k2 = 1:np;	
               amplitude(k1,k2)=amplitude(k1,k2)+c*PHold(k1, ip)*PPHold(ipp) ... 
                                *complex(cos(m*phi(k2)), -sin(m*phi(k2))) ...
                                * (1-2*mod((mabs-m)/2, 2));
            end;
         end;
      end;
   end;

   amplitude = conj(amplitude*complex(cos(fixphase(iPics)*pi/180.), sin(fixphase(iPics)*pi/180.)));    
                                     % put the matrix in the correct direction and fix the phase

   dist = abs(amplitude).^2;
   phase = atan2(imag(amplitude),real(amplitude));

   mx = max(max(abs(dist)));
   disp('mx'); disp(mx);
   titre=deblank(strjust(label(iPics,:),'left'));
   disp('Title');
   disp(titre);
   h=figure('visible','off');
   if (nvax > 1);
      vafac = 180/(nvax-1);
   else;
      vafac = 180.;
   end;

   for nva=1:nvax;
      va=(nva-1)*vafac;
      Construct_3D (Theta,Phi,dist,phase,mx,titre,va);
      filestem='';
      for ii=1:length(titre);
         if titre(ii) ~= ' ';
           filestem = strcat(filestem, titre(ii));
         end;
         if length(filestem) >= 15;
            break;
         end;
      end;
      if (nvax > 1);            
         filestem=strcat(filestem,'va',int2str(nva));
      end;
      disp('filestem');
      disp(filestem);
      saveas(h,strcat(filestem,'.eps'),'epsc');
      disp('Write files');
      disp(strcat(filestem,'.eps'));
   end;
end;


function yy=Construct_3D (thet,ph,rh,phase,m,tit,va)

% spherical coordinate equations
   r = rh.*sin(thet);
   x = r.*cos(ph);    % spherical coordinate equations
   y = r.*sin(ph);
   z = rh.*cos(thet);

% do the plot
   sh = surf(x,y,z,phase);
   set(sh,'AmbientStrength', 0.8);
% colormap(white);
   colormap(hsv);
   caxis([-pi pi]);
   colorbar;
   lh=light
   lightangle(-140,80);
%lighting phong
axis('square')

%par daut d'axis mise  l'helle automatique
   axis([-m m -m m -m m]);
   vaa=va-140;
   axis('on')
   view(vaa,30);

   xlabel('X','fontsize',14);ylabel('Y','fontsize',14);zlabel('Z','fontsize',14);
   title(tit,'fontsize',12);
%set(gca,'zTickLabel','-0.4|plan|0.1|0.2|0.4');
   set(gca,'fontsize',10);
