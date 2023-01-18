%% Assess floeDyn output file

%% Clean workspace
clear;
clc;
% close all;

%% Update path

addpath('../Process_scripts/');
addpath('../Analysis/');

%% Read in data

% outName = 'out_30p_50x50km_1237f_geo.h5'; 
% outPath = 'long_runs/';

outName = 'out_AKbB0.h5'; % long run no geo
outName = 'out_0UMab.h5'; % small floes
% outName = 'out_jVYVl.h5'; % 80% ice concentration run
outPath = '';


% outName = 'out_20p_50x50km_210f.h5';
% outName = 'out_10p_50x50km_99f.h5';
% outPath = 'sic_batches/FSD_300_7500_2/';

fo = readFloeOut( outName, outPath );
L = length(fo.state.time);
disp( fo.state.time(end)/86400 ) 

%% Plot floe trajectories

% sort by x location
[~,sInd] = sort( fo.state.xa(:,1) ); sInd = flipud(sInd);
% sInd = 1:fo.floes.numFloes;


ncol = max([fo.floes.numFloes,3]);
% ncol = L;
cols = colorcet( 'C7');
cols = cinterp(linspace(0,1,ncol),[0,1],cols);
% cols = colorcet( 'CBL2');
% cols = cinterp(linspace(0.15,0.95,ncol),[0,1],cols);
% cols = flipud(cols);

% cols = cols(end-ncol:end,:);

% cols = grey( linspace(0.15,0.85,ncol+1) );

psc = 1/1000; % plot scale (unit conversion m->km)
fontName = 'Open Sans';
fontSize = 14;

fH = figure(1); %clf;
hold on;
fH.Units = 'inches';
fH.Position([3,4]) = 6;

ax = gca;
hold on;
colororder( cols );
% plot( psc*fo.state.xa(:,end), psc*fo.state.ya(:,end), 'ko' );
% plot( psc*fo.state.xa.', psc*fo.state.ya.', LineWidth=1.5 );
plot( psc*fo.state.xp(sInd,:).', psc*fo.state.yp(sInd,:).', '.' );

xlim( psc*fo.chunk.win(1:2) );
ylim( psc*fo.chunk.win(3:4) );
daspect([1,1,1]);
xlabel('[km]')
ylabel('[km]')

ax.FontName = fontName;
ax.FontSize = fontSize;

drawnow;

%% Compare with forcing
% less of a "quick look"
%{
if ~isfield(fo,'ocean')
    forceName = '_input_forcing_SMgeo_5day.mat';
    forcePath = 'forcings/';
%     forceName = '_input_forcing.mat';
%     forcePath = '';    
    forcing = readForcing(forceName,forcePath);
    fo = floeflow( fo, forcing );
end

Uo = fliplr( fo.ocean.uo(:,2:end) + 1i*fo.ocean.vo(:,2:end)  ); 
Ui = fliplr( fo.state.u(:,2:end) + 1i*fo.state.v(:,2:end)  );
Im = fliplr( fo.state.im(:,2:end)  );
o2o = 0.2*[0,1];
uoss = linspace( o2o(1),o2o(2) );
uiss = seaIceSteadyState( 5.0e-3/1.4e-4, uoss );

fH = figure(2); %clf;
fH.Units = 'inches';
fH.Position([3,4]) = 6;

ax = gca;
hold on;
% plot( real(Uo(:)), real(Ui(:)),'.'); o2o = 0.2*[-1,1];
% plot( imag(Uo(:)), imag(Ui(:)),'.');

scH = scatter( abs(Uo(:)), abs(Ui(:)),[],log10(Im(:)),'filled');
scH.MarkerFaceAlpha = 0.65;
scH.MarkerEdgeColor = grey(0.25);
scH.MarkerEdgeAlpha = (0.25);
plot( o2o, o2o, 'k' );
plot( uoss, abs(uiss),Color=grey(0.5) );
xlim(o2o);
ylim(o2o);
caxis([5,10]);
daspect([1,1,1]);
xlabel('Floe-averaged current speed [m/s]')
ylabel('Sea ice translation speed [m/s]')

ax.FontName = fontName;
ax.FontSize = fontSize;

cmap = colorcet( 'L17', 'N', 20 );
cmap = cmap(6:end,:);
cmap = flipud( colorcet( 'L8', 'N', 20 ) );


colormap(cmap)
cb = colorbar;
ylabel(cb,'log_{10}( total received impulse )')
%}