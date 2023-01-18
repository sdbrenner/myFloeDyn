

%% Clean workspace

clearvars -except forcing;
clc;
% close all;

%% Load data

forcePath = 'forcings/';
forceName = '_input_forcing_SMgeo_long.mat';
if ~exist('forcing','var')
    forcing = readForcing(forceName,forcePath);
end

%% Make spectrogram

nx = length(forcing.x);
ny = length(forcing.y);

indx = 1:nx;
indy = floor(nx/2):nx;
[N,M,K] = size(forcing.u);
clear Zr;
for k = 1:K
    U = forcing.u(:,:,k)+1i*forcing.v(:,:,k);
    [kr,Zr(:,k)] = shellSpectra2D( forcing.x(indx), forcing.y(indy), U(indy,indx) );
end

%%

fH = figure(1); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [13.33,3.25];

cols = colorcet('CBTL1','N',ceil(1.05*K/18) );
colororder(cols);
% colormap(cols);

subplot(1,3,1);

loglog( kr,Zr(:,1:18:end), LineWidth=2  );
grid on;
% cb = colorbar;
dLabels = sprintfc('t=%2g days',forcing.t(1:18:end)/86400);
legH = legend(dLabels,Location='southwest',NumColumns=2,FontSize=6);

% figure(2); clf;

ax = subplot(1,3,2:3);

cmap = colorcet('CBL2','N',64 );
colormap(ax,cmap);
pcolor( forcing.t/86400, kr, log10(Zr) )
% shading flat;
ax.YScale = 'log';
% ax.CLim = [-9,-5.5];
% ax.CLim = [];
cb = colorbar;
