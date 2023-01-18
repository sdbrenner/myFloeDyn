%% Compare floe-averaged and floe-center-of-mass ocean properties


%% Clean workspace
clearvars -except forcing;
clc;


%% Load forcing data

forcePath = 'forcings/';
forceName = '_input_forcing_SMgeo_5day.mat';
if ~exist('forcing','var')
    forcing = readForcing(forceName,forcePath);
end

%% Load floe data

fDir = 'sic_batches/FSD_300_7500_2_geo/';
files = dir(['../../Floe_Cpp/io/outputs/',fDir,'*.h5']);
N = length(files);
n = 3;
fName = strcat(fDir,files(n).name);
fO = readFloeOut( files(n).name,fDir );
floeOut = floeflow( fO, forcing );

%% Plot velocity correlations

[M,L] = size( floeOut.ocean.com.uo );
fSzRange = [min(floeOut.floes.floeArea),max(floeOut.floes.floeArea)];
fSize = repmat( floeOut.floes.floeArea.', [1,L] );
msSzRange = ([9,200]);
msSize = interp1( fSzRange,msSzRange, fSize );
VAR = ( floeOut.ocean.var.uo(:) + floeOut.ocean.var.vo(:) );

fH = figure(1); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [7,6];
ax = gca;
hold on;
o2o = 0.2*[-1,1];
% dU = abs( floeOut.state.u(:) - floeOut.ocean.com.uo(:) );
scatter( floeOut.ocean.avg.uo(:), floeOut.ocean.com.uo(:),msSize(:),floeOut.ocean.var.uo(:),'filled',...
         MarkerEdgeColor=grey(0.15), MarkerEdgeAlpha=0.25, MarkerFaceAlpha=0.5 );
plot(o2o,o2o,'k');
ax.XLim = o2o;
ax.YLim = o2o;
ax.CLim = [0,2e-3];
cb = colorbar;
daspect([1,1,1]);
grid on;
xlabel('floe-averaged meridional ocean velocity [m/s]')
ylabel('floe-CoM meridional ocean velocity [m/s]')
ylabel(cb,'Velocity variance over floe [m^2/s^2]')

fH = figure(2); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [6,4];
ax = gca;
hold on;
scatter( sqrt(fSize(:)),VAR,9,'filled' );
ax.XScale = 'log';
% ax.YScale = 'log';
grid on;
xlabel('Floe size [m]')
ylabel('Velocity variance over floe [m^2/s^2]')
% xline(1500,'k--')

fH = figure(3); clf;
colororder( cbrewer2('-Set2') )
fH.Units = 'inches';
fH.Position([3,4]) = 6;
ax = gca;
hold on;
Ui = floeOut.state.ui(:,2:end);
UoAvg = floeOut.ocean.avg.uo(:,2:end);
UoCom = floeOut.ocean.com.uo(:,2:end);
msSize2 = msSize(:,2:end);
scatter( Ui(:),UoCom(:),msSize2(:),'filled',...
         MarkerEdgeColor=grey(0.15), MarkerEdgeAlpha=0.25, MarkerFaceAlpha=0.75 );
scatter( Ui(:),UoAvg(:),msSize2(:),'filled',...
         MarkerEdgeColor=grey(0.15), MarkerEdgeAlpha=0.25, MarkerFaceAlpha=0.75 );
plot(o2o,o2o,'k');
ax.XLim = o2o;
ax.YLim = o2o;
daspect([1,1,1]);
grid on;
xlabel('meridional ice velocity [m/s]')
ylabel('meridional ocean velocity [m/s]')
legend({'Floe-CoM','Floe-Avg'},Location='northwest');