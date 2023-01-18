%% Clean workspace
clearvars -except forcing;
clc;

%% Add paths

addpath('../Process_scripts/');

%% Load forcing data

forcePath = 'forcings/';
forceName = '_input_forcing_SMgeo_5day.mat';
if ~exist('forcing','var')
    forcing = readForcing(forceName,forcePath,calcFields=true);
end

%% Load floe data

% Data from long simulation:
% fName = 'out_30p_50x50km_1237f_geo.h5';
% fDir = 'long_runs/';

% Data from bacth simulations
fDir = 'sic_batches/FSD_300_7500_2_geo/';
files = dir(['../../Floe_Cpp/io/outputs/',fDir,'*.h5']);
N = length(files);
n = 3;

% Read in data
fName =  files(n).name;
floeOut = readFloeOut( fName,fDir );
% floeOut = floeflow( floeOut, forcing );

% % Using only common timesteps between forcing and FloeDyn output
% tol = 120; % timestep tolerance [seconds]
% [LIA,~] = ismembertol( floeOut.state.time, forcing.t,tol,'DataScale',1 );
% floeOut = subsetFloeOut( floeOut, LIA );

% Using only last timestep
% [~,L] = size( floeOut.state.ui );
% floeOut = subsetFloeOut( floeOut, L );

% Calculate 
floeOut = floeflow( floeOut, forcing, outputFloeMasks=true, outputDomainMask=true );
floeOut = makeDomainPoly(floeOut);
[floeOut,forcing] = calcIOBLDiss(floeOut,forcing);











%%
%{
%% Get ice velocity field

[M,L] = size(floeOut.floes.floeMasks);
[X,Y] = meshgrid( forcing.x, forcing.y );
W = range( floeOut.chunk.win(1:2));
H = range( floeOut.chunk.win(3:4));


Ui = NaN( [size(X),L] );
disp('Calculating full ice velocity field...')
% Loop through time
for l = 1:L
    UiL = NaN( size(X) );
    % Loop through floes
    for m = 1:M
        % translational ice velocity
        UiT = floeOut.state.ui(m,l) + 1i*floeOut.state.vi(m,l);
        % rotational ice velocity
        xc = floeOut.state.xp(m,l);
        yc = floeOut.state.yp(m,l);
        omega = floeOut.state.om(m,l);
        % Solve in periodic polar coords
%         if isPeriodicFloe(floeOut,m,l)
%             [TH,R] = cart2polPeriodic(X-xc,Y-yc,floeOut.chunk.win);
%         else
%             [TH,R] = cart2pol(X-xc,Y-yc);
%         end
%         UiR = omega.*R.*( -sin(TH) + 1i*cos(TH) );
        % Solve in cartesian coords (hopefully faster!)
        Xr = (X-xc) + W*((xc-X)>0.5*W) -W*((X-xc)>0.5*W);
        Yr = (Y-yc) + H*((yc-Y)>0.5*H) -H*((Y-yc)>0.5*H);
        UiR = omega*( -Yr + 1i*Xr );
        % total ice velocity
        iceMaskML = floeOut.floes.floeMasks{m,l};
        UiML = UiT + UiR;
        UiL(iceMaskML) = UiML(iceMaskML);       
    end
    Ui(:,:,l) = UiL;
    fprintf('timestep %2g of %2g...\n',l,L);
end


%% Calculate relative velocity, stress, and power/dissipation

fS = load('../MAT_Files/SMgeo_5day_iceVel.mat','ui','vi');
UiC = fS.ui + 1i*fS.vi;
UiC = UiC(:,:,1:L);

rho = 1025;
Cio = 5e-3;
f = 1.4e-4;
d = 1;

Uo = forcing.u(:,:,1:L) + 1i*forcing.v(:,:,1:L);
Ug = forcing.ug(:,:,1:L) + 1i*forcing.vg(:,:,1:L);
Urel = Uo-Ui;
Tio = rho*Cio*Urel.*abs(Urel);
diss = real(Tio.*conj(Urel)); 

UiSS = seaIceSteadyState( Cio/(f*d), Uo, Ug );
immoblDiss = rho*Cio*abs(Uo).^3;
steadyDiss = rho*Cio*abs(UiSS-Uo).^3;
cntuumDiss = rho*Cio*abs(UiC-Uo).^3;
excessDiss = diss./immoblDiss;


%% Average over floes or time

for l  = 1:L
    for m = 1:M
        iceMaskML = floeOut.floes.floeMasks{m,l};
        dissL = diss(:,:,l);
        steadyDissL = steadyDiss(:,:,l);
        floeOut.ocean.avg.diss(m,l) = mean( dissL(iceMaskML) );
%         floeOut.ocean.avg.excessDiss(m,l) = mean( dissL(iceMaskML)./steadyDissL(iceMaskML) );
        floeOut.ocean.avg.steadyDiss(m,l) = mean( steadyDissL(iceMaskML) );
    end
end


%}
%% Visualize
[M,L] = size(floeOut.state.ui);
[X,Y] = meshgrid( forcing.x, forcing.y );
l = L;

dissL = forcing.diss(:,:,l);
icePolyL = floeOut.chunk.icePoly{l};
icePolyL.Vertices = icePolyL.Vertices/1000;
imFrame = cell(1,L);

% cmap = colorcet('L8','N',20);
% cmap = colorcet('CBTL1','N',20);
cmap = flipud(colorcet('CBTL1','N',30)); cmap = cmap(1:20,:);

fH = figure(1); clf;
fH.Units = 'inches';
% fH.Position([3,4]) = 2*[5,4];
fH.Position([3,4]) = [6,5];
fH.Position([3,4]) = [3.25,3.1];
fH.Color = 'w';

SS = forcing.steadyDiss(:,:,l);
SS(~floeOut.chunk.iceMask{l}) = NaN;

ED = forcing.diss(:,:,l) - floeOut.chunk.iceMask{l}.*forcing.steadyDiss(:,:,l);
ED(ED<=0)=NaN;

ax = gca;
hold on;
% pH1 = pcolor( X/1000,Y/1000, log10(SS) ); 
pH2 = pcolor( X/1000,Y/1000, log10(forcing.diss(:,:,l)) );
% pH3 = pcolor( X/1000,Y/1000,log10(ED) ); 
% pgH = plot( icePolyL,FaceColor='none',EdgeColor=grey(0.5) );
% [~,cH] = contour(X/1000,Y/1000,abs(Uo(:,:,l)),...
%                  0.1*[1,1], Color=grey(0.75) );
ax.XLim = forcing.x([1,end])/1000;
ax.YLim = forcing.y([1,end])/1000;
ax.CLim = 0.003*[0,1];
ax.CLim = [-6,-1];
ax.FontSize = 8;
% ax.Color = 'k';

ax.XTickLabel = ''; ax.YTickLabel = '';
ax.Position = [0.025,0.025,0.95,0.95];

daspect([1,1,1]);
colormap( cmap );
% xlabel('[km]');
% ylabel('[km]');
% cb = colorbar;%('northoutside');
% cb.FontSize = 16;
% ylabel(cb,'log_{10}( \int\rho\epsilon{}dz )   [W/m^2]') 


%%
% for l = L
%     icePolyL = floeOut.chunk.icePoly{l};
%     icePolyL.Vertices = icePolyL.Vertices/1000;
%     pH1.CData = log10( abs(forcing.steadyDiss(:,:,l)) );
%     pH2.CData = log10( abs(forcing.diss(:,:,l)) );
%     pgH.Shape = icePolyL;
% %     cH.ZData = abs(Uo(:,:,l));
%     drawnow;
%     imFrame{l} = getframe(gcf);
% end

% Save gif
% imFrame(1) = [];
% gifPath = '~/Documents/Brown/Progress_meetings/Figures/2023/';
% gifFname = strcat(gifPath,fName(1:end-3),'notilt.gif');
% make_gif(gifFname,imFrame);


%%
%{

fH = figure(2); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [18,8];
fH.Color = 'w';

for l = L

icePolyL = floeOut.chunk.icePoly{l};
icePolyL.Vertices = icePolyL.Vertices/1000;

ax = subplot(2,3,1); cla;
hold on;
pH1 = pcolor( X/1000,Y/1000, log10(immoblDiss(:,:,l)) ); 
pgH1 = plot( icePolyL,FaceColor='none',EdgeColor=grey(0.5) );
shading flat;
ax.XLim = forcing.x([1,end])/1000;
ax.YLim = forcing.y([1,end])/1000;
ax.CLim = [-6,-1]; colormap( ax,cmap );
ax.FontSize = 12;
cb = colorbar;
daspect([1,1,1]);
xlabel('[km]');
ylabel('[km]');

ax = subplot(1,3,2); cla;
hold on;
pH2 = pcolor( X/1000,Y/1000, log10(diss(:,:,l)) ); 
pgH2 = plot( icePolyL,FaceColor='none',EdgeColor=grey(0.5) );
shading flat;
ax.XLim = forcing.x([1,end])/1000;
ax.YLim = forcing.y([1,end])/1000;
ax.CLim = [-6,-1]; colormap( ax,cmap );
ax.FontSize = 12;
cb = colorbar;
daspect([1,1,1]);
xlabel('[km]');
ylabel('[km]');

ax = subplot(2,3,3); cla;
hold on;
pH3 = pcolor( X/1000,Y/1000, real(log10(diss(:,:,l)./immoblDiss(:,:,l))) ); 
pgH3 = plot( icePolyL,FaceColor='none',EdgeColor=grey(0.5) );
shading flat;
ax.XLim = forcing.x([1,end])/1000;
ax.YLim = forcing.y([1,end])/1000;
ax.CLim = 4*[-1,1]; colormap( ax,cbrewer2('RdBu',21) );
% ax.CLim = [-6,-1]; colormap( ax,cmap );
ax.FontSize = 12;
cb = colorbar;
daspect([1,1,1]);
xlabel('[km]');
ylabel('[km]');



ax = subplot(2,3,4); cla;
hold on;
pH1 = pcolor( X/1000,Y/1000, log10(steadyDiss(:,:,l)) ); 
pgH1 = plot( icePolyL,FaceColor='none',EdgeColor=grey(0.5) );
shading flat;
ax.XLim = forcing.x([1,end])/1000;
ax.YLim = forcing.y([1,end])/1000;
ax.CLim = [-6,-1]; colormap( ax,cmap );
ax.FontSize = 12;
cb = colorbar;
daspect([1,1,1]);
xlabel('[km]');
ylabel('[km]');


ax = subplot(2,3,6); cla;
hold on;
pH3 = pcolor( X/1000,Y/1000, real(log10(diss(:,:,l)./steadyDiss(:,:,l))) ); 
pgH3 = plot( icePolyL,FaceColor='none',EdgeColor=grey(0.5) );
shading flat;
ax.XLim = forcing.x([1,end])/1000;
ax.YLim = forcing.y([1,end])/1000;
ax.CLim = 4*[-1,1]; colormap( ax,cbrewer2('RdBu',21) );
% ax.CLim = [-6,-1]; colormap( ax,cmap );
ax.FontSize = 12;
cb = colorbar;
daspect([1,1,1]);
xlabel('[km]');
ylabel('[km]');




drawnow;
end
linkaxes();


%%


bins = logspace( log10(300), log10(7500), 15 );
ind = 2:L;
dissML = floeOut.ocean.avg.excessDiss(:,ind);
floeSizeMat = floeOut.floes.floeSize.'*ones(1,size(dissML,2));
markerSizeMat = interp1( [300,7500],[2,12],floeSizeMat ).^2;
markerSize = 3.^2;
[~,binDiss25] = binFun( floeSizeMat(:), dissML(:), bins, binFunction=@(x) percentile(x,0.25) );
[~,binDiss50] = binFun( floeSizeMat(:), dissML(:), bins, binFunction=@(x) percentile(x,0.50) );
[binSz,binDiss75] = binFun( floeSizeMat(:), dissML(:), bins, binFunction=@(x) percentile(x,0.75) );


fH = figure(5); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [8,4];
cols = lines(8);

ax = gca;
hold on;
scatter( floeSizeMat(:), dissML(:),markerSize,...
         'filled',MarkerFaceAlpha=0.35, MarkerFaceCol=cols(1,:),...
                  MarkerEdgeAlpha=0.10, MarkerEdgeColor=grey(0.5) );
plot( binSz, binDiss50,'d',...
      MarkerSize=8,LineWidth=2,Color=colourShift(cols(1,:),-0.35) )
% errorbar( binSz, binDiss50, binDiss50-binDiss25, binDiss75-binDiss50,...
%           'd',MarkerSize=8,LineWidth=2,Color=colourShift(cols(1,:),-0.35) )
ax.XScale = 'log'; ax.YScale = 'log';
% ax.YLim = [1e-8,1e-2];
xline(1400,'k--',LineWidth=1 )

grid on;

xlabel('Floe size [m]');
ylabel('Excess dissipation ratio')
% ylabel('\int\epsilon{}dz   [W/m^2]') 


%%

ind = 2:L;
uoVar = floeOut.ocean.var.uo(:,ind) + floeOut.ocean.var.vo(:,ind);
dissML = floeOut.ocean.avg.excessDiss(:,ind);
floeSizeMat = floeOut.floes.floeSize.'*ones(1,size(dissML,2));
markerSizeMat = interp1( [300,7500],[2,12],floeSizeMat ).^2;


% o2o = [1e-10,1e-2];
% figure(8); clf; 
% ax = gca;
% hold on;
% scatter( uoVar(:), dissML(:),markerSizeMat(:),'filled',...
%          MarkerFaceAlpha=0.25, MarkerFaceCol=cols(1,:),...
%          MarkerEdgeAlpha=0.10, MarkerEdgeColor=grey(0.5) );
% plot( o2o, o2o,'-',Color=grey(0.5) );
% plot( o2o, rho*Cio*(o2o).^1.5,'k-.')
% ax.XLim = o2o; ax.YLim = o2o;
% ax.XScale = 'log'; ax.YScale = 'log';
% ax.FontSize = 12;
% % daspect([1,1,1]);
% grid on;
% xlabel('var( {\itu_o} ) [m^2/s^2]');
% ylabel('\langle\int\epsilon{}dz\rangle   [W/m^2]') 

%%
% 
% SIC(n) = floeOut.chunk.sic;
% meanDis(n) = mean( log10(diss(:)),'omitnan');


% function yn = isPeriodicFloe(floeOut,m,l)
%     vnX = floeOut.floes.floeXYT{m,l}(1,:).' ;
%     vnY = floeOut.floes.floeXYT{m,l}(2,:).' ;
%     win = floeOut.chunk.win;
%     yn = any( vnX<=win(1) | vnX>=win(2) | vnY<=win(3) | vnY>=win(4) );
% end
%}


%% Domain average?

nx = length(forcing.x); ny = length(forcing.y);
Uo = reshape( forcing.u(:,:,1:L)+1i*forcing.v(:,:,1:L), nx*ny,L);
Ui = reshape( forcing.ui+1i*forcing.vi, nx*ny,L);

UoDM = squeeze(mean(Uo));
UiDM = squeeze(mean(Ui,'omitnan')); UiDM(1) = NaN;

dissDM = 1025*5.5e-3*abs( UoDM-UiDM ).^3;

fH = figure(1); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [15,7];
ax = gca;
hold on;
qt = 0.003*(1:19);
quiver( qt,0*UoDM,real(UoDM),imag(UoDM),0,ShowArrowHead="off");
quiver( qt,0*UoDM,real(UiDM),imag(UiDM),0,ShowArrowHead="off");
daspect([1,1,1]);
