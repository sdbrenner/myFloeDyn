%% Investigate sea ice concentraiton changes with resolution

%% Clean workspace
clear;
clc;
close all;

%% Add path

addpath('../Process_scripts/');

%% Read in floe data

outPath = 'sic_batches/FSD_300_7500_2/';
outName = 'out_40p_50x50km_433f.h5';

fo = readFloeOut( outName, outPath );
L = length(fo.state.time);
fol = subsetFloeOut( fo, L );

%% Calculate binary mask for domain

domainX = 0:50:50e3;
domainY = 0:50:50e3;

fol  = makeDomainPoly(fol);
fol = makeFloeMask( fol, domainX, domainY );
[M,L] = size( fol.floes.floeMasks );

% Add domain floe mask (this is basically "free" once the floeMasks are calculated)
nx = length(domainX);
ny = length(domainY);
for l = 1:L
    floeMasksL = fol.floes.floeMasks(:,l);
    fml = reshape( full( [floeMasksL{:}] ), nx,ny,M);
    dml = sparse(any(fml,3));
    fol.chunk.iceMask{l} = dml;
end

iceMask = fol.chunk.iceMask{end};
imagesc(domainX,domainY,iceMask); shading flat;

%% Consider different ice concentrations

for binSize = [1,5,10,25,50,100,200] % [# of pts]

[y,x] = find(iceMask);
binEdges = 0:binSize:1000;
binX = 50*binCenters(binEdges);
N = histcounts2(x,y,binEdges,binEdges ).';
SIC = N/(binSize.^2);

figure(); clf;
ax = gca;
hold on;
imagesc(binX,binX,SIC); shading flat;
% pgH = plot( fol.chunk.icePoly{1} );
% pgH.FaceColor = 'none';
% pgH.LineWidth = 2;
daspect([1,1,1]);
colorbar;
ax.XLim = [0,50e3];
ax.YLim = [0,50e3];
ax.CLim = [0,1];

% colormap( colorcet('CBD2','N',10) );
colormap( cbrewer2('Blues',10) );
mean(SIC(:))
% waitforbuttonpress;

end




% Something else to try?
% A=[1;1;2;3;3;3];
% B=rand(6,10);
% [xx, yy] = ndgrid(A,1:size(B,2));
% C=accumarray([xx(:) yy(:)],B(:));