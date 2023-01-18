% Batch evaluate

%% Clean workspace
clear;
clc;
% close all;

%% Add 'process' path

addpath('../Process_scripts/');

%% Load data

batchCase = 'FSD_300_2500_2_geo';
fName = strcat('floeOut_',batchCase );
load( strcat('../MAT_Files/',fName) );
floeOutOrignal = floeOut;
%% Calculate correlations

N = length(floeOut);
R = NaN( 2,N );
SIC = NaN( 1,N );

for n = 1:N
    SIC(n) = floeOut(n).chunk.sic;
    Ui = floeOut(n).state.ui  + 1i*floeOut(n).state.vi;
    Uo = floeOut(n).ocean.avg.uo + 1i*floeOut(n).ocean.avg.vo;
%     Uiss = steadyUiInt( Uo,1.4e-4,5e-3 );

    rr = corrcoef( abs(Uo(:,2:end)), abs(Ui(:,2:end)) );
    R(1,n) = rr(1,2);
%     rr = corrcoef( real(Uiss(:,2:end)), real(Ui(:,2:end)) );
%     R(2,n) = rr(1,2);
%     rr = corrcoef( imag(Uiss(:,2:end)), imag(Ui(:,2:end)) );
%     R(3,n) = rr(1,2);

    rr = corrcoef( floeOut(n).state.om(:,2:end), floeOut(n).ocean.avg.zeta(:,2:end) );
    R(4,n) = rr(1,2);
end

%% Plot correlations

% fontName = 'Century Gothic';
fontName = 'Montserrat';
fontSize = 14;
mea = 1;
mfa = 0.5; 
ms = 10;
cols = lines(7);

fH = figure(1); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [5.5,5.5];

ax = subplot(2,1,1);
hold on;
scatter( SIC, R(1,:),ms^2,'filled',...
         MarkerEdgeColor=grey(0.5),MarkerEdgeAlpha=mea, MarkerFaceAlpha=mfa );
% scatter( SIC, R(2,:),(0.65*ms)^2,'filled',...
%          MarkerEdgeColor=grey(0.5),MarkerEdgeAlpha=mea, MarkerFaceAlpha=mfa );
% scatter( SIC, R(3,:),(0.65*ms)^2,'filled',...
%          MarkerEdgeColor=grey(0.5),MarkerEdgeAlpha=mea, MarkerFaceAlpha=mfa );
ax.XLim = [0,1];
ax.YLim = [0,1];
ax.YTick = 0:0.25:1;
ax.FontName = fontName;
ax.FontSize = fontSize;
grid on;
xlabel('Sea ice concentration')
% ylabel('Velocity correlation')
ylabel('Corr[ {\bfU}_i, \langle{\bfu}_o\rangle ]')

ax = subplot(2,1,2);
hold on;
scatter( SIC, R(4,:),ms^2,'filled',...
         MarkerEdgeColor=grey(0.5),MarkerEdgeAlpha=mea, MarkerFaceAlpha=mfa );
ax.XLim = [0,1];
ax.YLim = [0,1];
ax.YTick = 0:0.25:1;
ax.FontName = fontName;
ax.FontSize = fontSize;
grid on;
xlabel('Sea ice concentration')
% ylabel('Rotation correlation')
ylabel('Corr[ \omega_i, \langle\zeta_o\rangle ]')


% close all;

%% Plot

fH = figure(2); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [7.5,5.5];

floeOut = floeOutOrignal;
% floeOut = floeOutOrignal(1:2:end);
N = length(floeOut);


% W = 0.09;
W = 0.18;
H = fH.Position(3)/fH.Position(4)*W;
sx = (1-N*W)/(N+1);
sy = (1-3*H)/4;
fontSize = 8;

h = 1;
f = 1.4e-4;
C = 5e-3;
uos = linspace(0,0.2);
uis = seaIceSteadyState( C/(f*h), uos );

for n = 1:N
    % Extract plot data
    SIC(n) = floeOut(n).chunk.sic;
    [M,L] = size(floeOut(n).state.ui);
    plotInds = 2:L;
    imNorm = fliplr( floeOut(n).state.im(:,plotInds)./floeOut(n).floes.floeArea.' );
    Ui =     fliplr( floeOut(n).state.ui(:,plotInds)  + 1i*floeOut(n).state.vi(:,plotInds) );
    Uo =     fliplr( floeOut(n).ocean.avg.uo(:,plotInds) + 1i*floeOut(n).ocean.avg.vo(:,plotInds) );
    Zeta =   fliplr( floeOut(n).ocean.avg.zeta(:,plotInds) );
    Omega =  fliplr( floeOut(n).state.om(:,plotInds) );
    % Generate scatterplot marker sizes based on floe sizes
    msSizeRange = [2,10];
    msSize = interp1( floeOut(n).floes.floeSize([end,1]),msSizeRange,...
                      floeOut(n).floes.floeSize);
    szd = repmat(msSize.',1,size(Ui,2) );
    msSize = fliplr(szd);
    % Generate polyshape of final sea ice configuration
    fol = subsetFloeOut( floeOut(n), L );
    fol = makeDomainPoly( fol );

    % Plot: translation
    figure(2)
    ax = axes();
    o2o = 0.2*[0,1];
    hold on;
    scatter( abs(Ui(:)), abs(Uo(:)),szd(:).^2,log10(imNorm(:)),...
             'filled',MarkerFaceAlpha=0.5 );
    plot( o2o, o2o,'k-')
    plot( abs(uis), abs(uos),Color=grey(0.5) );
    ax.XLim = o2o;
    ax.YLim = o2o;
    ax.YTick = 0:0.1:0.2; ax.XTick = 0:0.1:0.2;
    ax.CLim = [-3,5];
    X = 1.25*sx+(n-1)*(W+sx);
    Y = sy+(H+sy);
    ax.Position = [X,Y,W,H];
    ax.FontName = fontName;
    ax.FontSize = fontSize;
    grid on; grid minor;
    daspect([1,1,1]);
    xlabel('||{\bfU}_i||   [m/s]')
    if n==1, ylabel('||\langle{\bfu}_o\rangle||   [m/s]'); end

    % Plot: rotation
    ax = axes();
    o2o = [-1.5,2.5]; %2*[-1,1];
    hold on;
    scatter( Omega(:)/f, Zeta(:)/f, szd(:).^2, log10(imNorm(:) ) ,...
             'filled',MarkerFaceAlpha=0.5 );
    plot( o2o, 2*o2o,'k-',0,0,'k+',MarkerSize=10);
    ax.XLim = o2o;
    ax.YLim = o2o; 
    ax.YTick = -3:1:3; ax.XTick = -3:1:3;
    ax.CLim = [-3,5];
    X = 1.25*sx+(n-1)*(W+sx);
    Y = sy;
    ax.Position = [X,Y,W,H];  
    ax.FontName = fontName;
    ax.FontSize = fontSize;    
    grid on;
    daspect([1,1,1]);
    xlabel('\omega_i/f')
    if n==1, ylabel('\langle\zeta_o\rangle/f'); end

    % Plot: ice assembly
    ax = axes();
    hold on;
%     ax.Color = [0.3,0.4,0.6];
    pgH = plot( fol.chunk.icePoly{end} );
%     pgH.FaceAlpha = 1; pgH.FaceColor=grey(0.7); pgH.EdgeColor='none';
    o2o = 50e3*[0,1];
    ax.XLim = fol.chunk.win(1:2);
    ax.YLim = fol.chunk.win(3:4); 
    ax.XTick = []; 
    ax.YTick = [];
    X = 1.25*sx+(n-1)*(W+sx);
    Y = sy+2*(H+sy);
    ax.Position = [X,Y,W,H];  
    ax.FontName = fontName;
    ax.FontSize = fontSize;    
    daspect([1,1,1]);
    titLab = sprintf('SIC = %2.0d%%',round(100*SIC(n)) );
    title( titLab );

    cmap = colorcet('L17','N',20);
    colormap( cmap(5:end,:) );

    drawnow;
end





