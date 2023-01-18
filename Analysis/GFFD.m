%% GFFD
% Calculate the floe-floe-distance-distribution function (Herman 2011)

%% Clean workspace
clear;
clc;
% close all;


%% Add paths

addpath('../Process_scripts/');

%% Load floe data

% Data from bacth simulations
fDir = 'sic_batches/FSD_300_7500_2/';
files = dir(['../../Floe_Cpp/io/outputs/',fDir,'*.h5']);
N = length(files);
n = 5;

% Read in data
fName =  files(n).name;
floeOut = readFloeOut( fName,fDir );
floeOut = calcFloeXYT( floeOut );

% Use only first and last timesteps
[~,L] = size( floeOut.state.ui );
floeOut = subsetFloeOut( floeOut, [1,L] );
floeOut = makeDomainPoly( floeOut );

%% Calculate minimum floe-floe distances and GFFD

% Extract data from floeOut
[M,L] = size( floeOut.floes.floeXYT );
W = range(floeOut.chunk.win(1:2));
H = range(floeOut.chunk.win(3:4));

% PDF bin properties:
minDist = 0; % semi-arbitrary choice 
% maxDist = 0.5*sqrt( W^2 + H^2 );
maxDist = min([W,H])/2; % for maxDist larger than this, need to account for periodic area in normalization
numBins = 100;

% Define bin edges for PDF: log-spaced
% binEdges = logspace( log10(minDist), log10(maxDist), numBins+1 );
% bins = 10.^binCenters(log10(binEdges));

% Define bin edges for PDF: linear-spaced
binEdges = linspace( minDist,maxDist, numBins+1 );
bins = binCenters( binEdges );


% Calculate normalization factor
% (following Luding & Strauß, 2001; eqn. 4)
binSizes = diff(binEdges);
V = pi*maxDist.^2;
Vr = pi*(2*binEdges(1:end-1) + binSizes).*binSizes;
numPairs = M*(M-1)/2;
normFactor = (1/numPairs)*(V./Vr);

% Pre-allocate
dists = NaN( 0.5*M*(M-1), L );
gffd = NaN(numBins,L);

% Loop through time
for l = 1:L
    % Loop through floe-pairs
    dmin = NaN(M,M);
    for i = 1:M
        for j = 1:i-1 % only need to cover upper triangle of matrix
            % all distances pairs between two floes   
            fXYi = floeOut.floes.floeXYT{i,l};
            fXYj = floeOut.floes.floeXYT{j,l};
            d = toroidalDistIJ( fXYi, fXYj, floeOut.chunk.win); 
            % minimum distance
            [dmin(i,j), idx] = min(d(:));
%             [row, col] = ind2sub(size(d), idx);
%             figure(1); 
%             hold on;
%             plot( fXYi(1,:),fXYi(2,:),fXYj(1,:),fXYj(2,:) ); 
%             plot( [fXYi(1,col),fXYj(1,row)],[fXYi(2,col),fXYj(2,row)] ); 
%             xlim(floeOut.chunk.win([1,2]));
%             ylim(floeOut.chunk.win([3,4]));
%             daspect([1,1,1]);
%             drawnow;
%             waitforbuttonpress;
        end
    end
    dminvect = dmin(:);
    dminvect(isnan(dminvect)) = [];
    dists(:,l) = dminvect;
    distCounts = histcounts( dists(:,l),binEdges,'normalization','count' );
    normFactCorr = numPairs/sum(distCounts); % accounts for distance outside of bin ranges
    gffd(:,l) = distCounts .* normFactor*normFactCorr; % normalization 
%     gffd(:,l) = histcounts( dists(:,l),binEdges,'normalization','pdf' );
end

%%

fH = figure(1); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [6,6];

col = lines(8);

ax = subplot(2,2,1); cla;
icePoly1 = floeOut.chunk.icePoly{1};
icePoly1.Vertices = icePoly1.Vertices/1000;
plot( icePoly1,FaceColor=col(1,:));
ax.XLim =  floeOut.chunk.win(1:2)/1000 ;
ax.YLim =  floeOut.chunk.win(3:4)/1000 ;
daspect([1,1,1]);
title('t = 0')

ax = subplot(2,2,2); cla;
icePolyL = floeOut.chunk.icePoly{L};
icePolyL.Vertices = icePolyL.Vertices/1000;
plot( icePolyL,FaceColor=col(L,:));
ax.XLim =  floeOut.chunk.win(1:2)/1000 ;
ax.YLim =  floeOut.chunk.win(3:4)/1000 ;
daspect([1,1,1]);
title('t = 72 h')

ax = subplot(2,2,3:4);
loglog( bins, gffd, '.-',LineWidth=1.5 );
hold on;
% ax.YScale = 'linear';
xlabel('Floe-floe distance [m]')
ylabel('g_{ffd}(d)');
grid on;
ax.FontSize = 12;
% xline(W/2,'k');
yline( 1 ,'k--')

%% EMBEDDED FUNCTIONS

function dists = toroidalDistIJ(xyi,xyj,win)
    
    X = win(1); Y = win(3);
    W = range(win(1:2)); H = range(win(3:4));

    % transform from X+[0,W]->[0,W],Y+[0,H]->[0,H]
    xyi = xyi-[X;Y];
    xyj = xyj-[X;Y];
    
    % Calculate regular dx dy for each point set
    dx = abs(xyi(1,:)-xyj(1,:).');
    dy = abs(xyi(2,:)-xyj(2,:).');
    
    % Modify dx dy
    dx(dx>0.5*W) = W-dx(dx>0.5*W);
    dy(dy>0.5*H) = H-dy(dy>0.5*H);
    
    % now calculate distance
    dists = sqrt( dx.^2 + dy.^2 );

end
