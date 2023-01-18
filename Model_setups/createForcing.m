%% CREATE ANNA FORCING
% Load, visualize, and restructure model output data

%% Clean workspace and load data 
clear;
clc;
close all;

load('../MAT_Files/3D_Standard_highstrat_surface_velocity.mat');
load('../MAT_Files/3D_Standard_highstrat_surface_fields.mat','ETAN','RHO');

%% Tranpose matrices

% Correct SSH at first column
ETAN(:,1,:) = ETAN(:,2,:);

% Transpose (permute) U,V,SSH matrices
ETAN = permute(ETAN,[2 1 3]);
RHO = permute(RHO,[2 1 3]);
U = permute(U,[2 1 3]);
V = permute(V,[2 1 3]);




%% Calculate velocity fields

[N,M,K] = size(U);

f = 1.4e-4;
g = 9.81; 

Ug = NaN( N,M,K ); Vg = NaN( N,M,K );
for k = 1:K 
    [etax,etay] = gradient( ETAN(:,:,k),X,Y);   
%     ugck = g/(1i*f)*(etax+1i*etay);
    Ug(:,:,k) = -g/f*etay;
    Vg(:,:,k) = +g/f*etax;
end


% Speeds:
spd = abs(U+1i*V);
geospd = abs(Ug+1i*Vg);

%% Animate
%{
fH = figure(1); clf;
ax = gca;
hold on;

k = 1;
pH = pcolor( X,Y,spd(:,:,k)); shading interp;
qInd = 1:20:length(X);
qH = quiver( X(qInd),Y(qInd),U(qInd,qInd,k),V(qInd,qInd,k),1,'k');
daspect([1,1,1]);
ax.XLim = X([1,end]);
ax.YLim = Y([1,end]);
xline(25e3,'k',LineWidth=1.5); 

cmap = flipud(colorcet('L17','N',20));
% cmap = cbrewer2('RdBu',21);
colormap( cmap );
cb = colorbar;


for k = 1:K %(6*8)
    [zetak,~] = curl( X,Y,V(:,:,k),U(:,:,k) );
    divk = divergence( X,Y,V(:,:,k),U(:,:,k) );
    pH.CData = geospd(:,:,k); caxis(0.2*[0,1]);
%     pH.CData = ETAN(:,:,k); caxis(0.01*[-1,1]);
%     pH.CData = spd(:,:,k); caxis(0.2*[0,1]);
%     pH.CData = zetak; caxis(2*f*[-1,1]);
    qH.UData = -V(qInd,qInd,k); qH.VData = -U(qInd,qInd,k);
    drawnow;
end
%}
%% Compare geostrophic and ageostrophic velocity


fH = figure(9); clf;
fH.PaperUnits = 'inches';
fH.PaperPosition([3,4]) = [13.33,7.5];
fH.Color = 'w';
dpi = 150;
fH.Position([3,4]) = fH.PaperPosition([3,4])*dpi;

qInd = 1:20:length(X);
Xkm = X/1000; Ykm = Y/1000; 
asc = 5; % "arrow scale factor"
CL = 0.25*[0,1];
k = 48;
cmap = flipud(colorcet('L17','N',20));

ax = subplot(2,3,1);
hold on;
pH1 = pcolor( Xkm,Ykm,spd(:,:,k)); shading interp;
qH1 = quiver( Xkm(qInd),Ykm(qInd), ...
              asc*U(qInd,qInd,k), asc*V(qInd,qInd,k), ...
              'k','AutoScale','off');
daspect([1,1,1]);
ax.XLim = Xkm([1,end]);
ax.YLim = Ykm([1,end]);
ax.CLim = CL;
% xline(25e3,'k',LineWidth=1.5); 
cb = colorbar;
colormap( ax,cmap );
title('Total velocity');


ax = subplot(2,3,2);
hold on;
pH2 = pcolor( Xkm,Ykm,geospd(:,:,k)); shading interp;
qH2 = quiver( Xkm(qInd),Ykm(qInd), ...
              asc*Ug(qInd,qInd,k),asc*Vg(qInd,qInd,k), ...
              'k','AutoScale','off');
daspect([1,1,1]);
ax.XLim = Xkm([1,end]);
ax.YLim = Ykm([1,end]);
ax.CLim = CL;
% xline(25e3,'k',LineWidth=1.5); 
cb = colorbar;
colormap( ax,cmap );
title('"Geostrophic" velocity');


ax = subplot(2,3,3);
Ua = U-Ug; Va = V-Vg; 
ageospd = abs( Ua+1i*Va );
hold on;
pH3 = pcolor( Xkm,Ykm,ageospd(:,:,k)); shading interp;
qH3 = quiver( Xkm(qInd),Ykm(qInd), ...
              asc*Ua(qInd,qInd,k),asc*Va(qInd,qInd,k), ...
              'k','AutoScale','off');
daspect([1,1,1]);
ax.XLim = Xkm([1,end]);
ax.YLim = Ykm([1,end]);
ax.CLim = CL;
% xline(25e3,'k',LineWidth=1.5); 
cb = colorbar;
colormap( ax,cmap );
title('"Ageostrophic" (residual) velocity');

linkaxes();

%%
cmap = cbrewer2('RdBu',21);
[zetak,~] = curl( X,Y,U(:,:,k),V(:,:,k) );
divk = divergence( X,Y,U(:,:,k),V(:,:,k) );

ax = subplot(2,3,4);
hold on;
pH4 = pcolor( Xkm,Ykm,divk); shading interp;
daspect([1,1,1]);
ax.XLim = Xkm([1,end]);
ax.YLim = Ykm([1,end]);
ax.CLim = 3d-4*[-1,1];
% xline(25e3,'k',LineWidth=1.5); 
cb = colorbar;
colormap( ax,cmap );
title('Divergence');


ax = subplot(2,3,5);
hold on;
pH5 = pcolor( Xkm,Ykm,zetak/f); shading interp;
daspect([1,1,1]);
ax.XLim = Xkm([1,end]);
ax.YLim = Ykm([1,end]);
ax.CLim = 3*[-1,1];
% xline(25e3,'k',LineWidth=1.5); 
cb = colorbar;
colormap( ax,cmap );
title('Vorticity');


ax = subplot(2,3,6);
hold on;
pH6 = pcolor( Xkm,Ykm,ETAN(:,:,k)); shading interp;
daspect([1,1,1]);
ax.XLim = Xkm([1,end]);
ax.YLim = Ykm([1,end]);
ax.CLim = 0.01*[-1,1];
% xline(25e3,'k',LineWidth=1.5); 
cb = colorbar;
colormap( ax,cmap );
title('SSH');



%% Animate 
asc = 5;
imFrame = cell(1,K);

for k = 1:K
    % Calculate vorticity and divergence
    [zetak,~] = curl( X,Y,U(:,:,k),V(:,:,k) );
    divk = divergence( X,Y,U(:,:,k),V(:,:,k) );

    % Update shadings: speeds
    pH1.CData = spd(:,:,k);
    pH2.CData = geospd(:,:,k);
    pH3.CData = ageospd(:,:,k);
    % Update vectors
    qH1.UData = asc*U(qInd,qInd,k); qH1.VData = asc*V(qInd,qInd,k);
    qH2.UData = asc*Ug(qInd,qInd,k); qH2.VData = asc*Vg(qInd,qInd,k);
    qH3.UData = asc*Ua(qInd,qInd,k); qH3.VData = asc*Va(qInd,qInd,k);
    % Update shadings: other
    pH4.CData = divk;
    pH5.CData = zetak/f;
    pH6.CData = ETAN(:,:,k);

    % Update plot
    drawnow;
    imFrame{k} = getframe(fH);
end

make_gif('forcing.gif',imFrame);





%% Organize and save data output

forcing.x = X.';
forcing.y = Y.';
forcing.t = T.'*86400;
forcing.u = U; 
forcing.v = V;
forcing.ug = Ug; 
forcing.vg = Vg;

forceName = '_input_forcing_SMgeo_long.mat';
forcePath = '../../Floe_Cpp/io/inputs/forcings/';
fName = strcat(forcePath,forceName);
save(fName,'-struct','forcing');
fprintf('Saved %s\n',fName);


%% Create subset beginning on day 25 of the simulation

ind = T>=25;
forcing.t = forcing.t(ind);
forcing.t = forcing.t-forcing.t(1);
forcing.u = forcing.u(:,:,ind);
forcing.v = forcing.v(:,:,ind);
forcing.ug = forcing.ug(:,:,ind);
forcing.vg = forcing.vg(:,:,ind);

forceName = '_input_forcing_SMgeo_5day.mat';
forcePath = '../../Floe_Cpp/io/inputs/forcings/';
fName = strcat(forcePath,forceName);
save(fName,'-struct','forcing');
fprintf('Saved %s\n',fName);


%%

forcing.X = X.';
forcing.Y = Y.';
% forcing.T = T.'*86400;
forcing.U = U(:,:,end); 
forcing.V = V(:,:,end);
forcing.RHO = RHO(:,:,end);
forcing.ETAN = ETAN(:,:,end);

fName = 'fTest.mat';
save(fName,'-struct','forcing');
