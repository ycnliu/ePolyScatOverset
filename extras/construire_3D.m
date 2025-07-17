function yy=construire_3D (thet,ph,rh,m,tit, viewaz, viewel)

% spherical coordinate equations
r = rh.*sin(thet);
x = r.*cos(ph);    % spherical coordinate equations
y = r.*sin(ph);
z = rh.*cos(thet);

% do the plot
surf(x,y,z)
% colormap(white);
colormap(autumn);
%colorbar;
%light
%lighting phong
lightangle(viewaz,viewel);
axis('square');

%par daut d'axis mise  l'helle automatique
axis([-m m -m m -m m]);

axis('on');
%view(-140,30);
%view(-140,80);
view(viewaz,viewel);

xlabel('X','fontsize',28);ylabel('Y','fontsize',28);zlabel('Z','fontsize',28);
% title(tit,'fontsize',12);
%set(gca,'zTickLabel','-0.4|plan|0.1|0.2|0.4');
set(gca,'fontsize',20);
% pwh=get(gcf,'PaperSize');
% disp('width and height');
% disp(pwh(1));
% disp(pwh(2));
