%% Example
% Example of data processing workflow

% Clean workspace
clear;
clc;

% Forcing
forceName = '_input_forcing_SM_5day.mat';
forcePath = 'forcings/';
forcing = readForcing( forceName, forcePath );

%% Floes
outPath = 'sic_batches/FSD_300_7500_2/';
outName = 'out_30p_50x50km_322f.h5';

floeOut = readFloeOut( outName, outPath );
% floeOut = makeFloeMask( floeOut,forcing.x,forcing.y );
% floeOut = floeflow( floeOut, forcing );
tic
% floeOut = makeDomainPoly( floeOut );
runTime = toc;

%%
figure(1); clf;
[M,L] = size(floeOut.floes.floeXYT);

l = L;
floeMasksL = floeOut.floes.floeMasks(:,l);
nx = length(forcing.x);
ny = length(forcing.y);
A = reshape( full( [floeMasksL{:}] ), nx,ny,M);
dml = sparse(any(A,3));
imagesc( forcing.x,forcing.y,dml );

hold on;
for m = 1:M
    fXYT = floeOut.floes.floeXYT{m,l};
    plot( fXYT(1,:), fXYT(2,:), LineWidth = 2 );
end
daspect([1,1,1]);
xlim( floeOut.chunk.win(1:2) );
ylim( floeOut.chunk.win(3:4) );
%%
figure(2); clf;
hold on;
Ui = floeOut.state.u(:,l) + 1i*floeOut.state.v(:,l);
Uo = floeOut.ocean.uo(:,l) + 1i*floeOut.ocean.vo(:,l);
% plot( floeOut.state.u(:,l), floeOut.ocean.uo(:,l), '.')
% plot( floeOut.state.v(:,l), floeOut.ocean.vo(:,l), '.')
plot( abs(Ui),abs(Uo),'.' );
daspect([1,1,1]);