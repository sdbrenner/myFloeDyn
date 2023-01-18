
clearvars -except forcing;
clc;

forcePath = 'forcings/';
forceName = '_input_forcing_SMgeo_5day.mat';
if ~exist('forcing','var')
    forcing = readForcing(forceName,forcePath);
end

[N,M,K] = size( forcing.u );
k = K;
X = forcing.x;
Y = forcing.y;
uo = forcing.u(:,:,k) + 1i*forcing.v(:,:,k);
ug = forcing.ug(:,:,k) + 1i*forcing.vg(:,:,k);
[zeta,~] = curl(X,Y,real(uo),imag(uo) );

Cstar = 5e-3/(1.4e-4);
[ui1,err1] = seaIceSteadyState( Cstar,uo,0 );
[ui2,err2] = seaIceSteadyState( Cstar,uo,ug );
%%
div = @(U) divergence(X,Y,real(U),imag(U) );

fH = figure(1); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [13.33,4];
colormap( cbrewer2('RdBu',21) );

ax = subplot(1,3,1);
pcolor( X/1000,Y/1000,div(uo) );
ax.XLim = [0,50]; ax.YLim = [0,50];
ax.CLim = 2e-4*[-1,1];
ax.XTick = 0:10:50;
ax.YTick = 0:10:50;
daspect([1,1,1]);
ax.Position(2) = 0.22;

ax = subplot(1,3,2);
hold on;
pcolor( X/1000,Y/1000,div(ui1) );
ax.XLim = [0,50]; ax.YLim = [0,50];
ax.CLim = 2e-4*[-1,1];
ax.XTick = 0:10:50;
ax.YTick = 0:10:50;
daspect([1,1,1]);
ax.Position(2) = 0.22;

cb = colorbar('southoutside');
cb.Position([1,3]) = ax.Position([1,3]);
cb.Position([2,4]) = [0.11,0.035];
xlabel(cb,'divergence [1/s]')

ax = subplot(1,3,3);
pcolor( X/1000,Y/1000,div(ui2) );
ax.XLim = [0,50]; ax.YLim = [0,50];
ax.CLim = 2e-4*[-1,1];
ax.XTick = 0:10:50;
ax.YTick = 0:10:50;
daspect([1,1,1]);
ax.Position(2) = 0.22;



%%
figure(2); clf;

subplot(1,2,1);
pcolor(X,Y,err1);
caxis([0,1e-3]);
daspect([1,1,1]);
colorbar('southoutside');

subplot(1,2,2);
pcolor(X,Y,err2);
caxis([0,1e-3]);
daspect([1,1,1]);
colorbar('southoutside');
