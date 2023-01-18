%% SCALE DEPENDANCE OF IDEAL "FLOE-AVERAGED" FORCING

%% Clear workspace

clear;
clc;
close all;

%% Add paths

addpath('../Process_scripts/');

%% Load forcing data

% forceName = '_input_forcing_SMgeo_long.mat';
forceName = '_input_forcing_SMgeo_5day.mat';
forcePath = 'forcings/';
forcing = readForcing(forceName,forcePath);

% Floe-shape inventory
load('../../Floe_Cpp/io/inputs/Biblio_Floes.mat');

%% Extract data for given timestep

k = 1;
Uo = forcing.u(:,:,k) + 1i*forcing.v(:,:,k);
[zeta,~] = curl( forcing.x, forcing.y, real(Uo), imag(Uo) );


%% Design sampling distribution

minSize = 300;
maxSize = 30e3;
numSizes = 35;
numFloes = 10;

% logarithmic floe-size distribution
sf = logspace(log10(minSize),log10(maxSize),numSizes);
sf = repmat(sf,[numFloes,1]);


%% Randomly sample domain

% Window and boundary size
X = forcing.x(1); 
Y = forcing.y(1); 
W = range(forcing.x); 
H = range(forcing.y); 
dx = mean(diff(forcing.x)); dy = mean(diff(forcing.y));
position = [X,Y,W,H];
window = [X+W*[0,1],Y+H*[0,1]];
boundary = [ X+[0;0;1;1;0]*W , Y+[0;1;1;0;0]*H];

% Initialize "floeOut" structure 
fo.chunk.win = window;  
fo.floes.floeSize = sf;
fo.floes.floeArea = sf.^2;
fo.state.time = NaN(1,numSizes);
fo.state.PBC = 1;

for l = 1:numFloes
    for m = 1:numSizes
        % Generate random floe position and orientation
        floeResolution=150;
        floeK = floeInventory( G,sf(n),floeResolution );
        % Create equivilant "floeOut" structure
        fo.floes.floeShapes(m,l) = floeK;
        fo.state.xa(m,l) = W*rand(1)+X;
        fo.state.ya(m,l) = H*rand(1)+Y;
        fo.state.th = 2*pi*rand(1); 
        fo.state.xp(m,l) = fo.state.xa(m,l);
        fo.state.yp(m,l) = fo.state.ya(m,l);       
    end
end

fo = floeflow(fo,forcing);



%%
%{









% Pre-allocate
UFloe = NaN(size(sf));
zetaFloe = NaN(size(sf));

% Loop through floes
for n = 1:numel(sf)
        % Generate random floe position and orientation
        floeResolution=150;
        floeK = floeInventory( G,sf(n),floeResolution );
        theta = 2*pi*rand(1);
        fsXY(1) = W*rand(1)+X;
        fsXY(2) = H*rand(1)+Y;
        R = [ cos(theta), -sin(theta) ; 
              sin(theta),  cos(theta)  ];

        % Put floe into position
        floePos = (floeK*R)+fsXY;
        floeXY{1} = (floeK*R)';

        % Create equivilant "floeOut" structure
        thisFloe.state.PBC = 1;
        thisFloe.chunk.win = window;
        thisFloe.floes.floeXYT{1} = floePos.';

        % Get binary mask for floe
        thisFloe = makeFloeMask( thisFloe, forcing.x, forcing.y );
        floeMask = thisFloe.floes.floeMasks{1};

        % Average ocean velocity and vorticity under floe
        UFloe(n) = mean(Uo(floeMask));
        zetaFloe(n) = mean(zeta(floeMask));

        % Visualize (as a check)
%         figure(1); clf;
%         pcolor( floeMask ); 
%         caxis(0.2*[0,1]);
%         daspect([1,1,1]);
%         drawnow;
end



%% Plot output


fH = figure(1); clf;

ax = subplot(2,1,1);
hold on;
plot( sf(:), abs(UFloe(:)).^2,'.',Color=grey(0.8) );
plot( mean(sf), mean(abs(UFloe).^2),'k-',LineWidth=1.5 );
ax.XScale = 'log';
ax.YScale = 'log';
ax.XLim = [50,50e3];
ax.YLim = [0.01,0.25].^2;


ax = subplot(2,1,2);
hold on;
plot( sf(:), abs(zetaFloe(:)).^2,'.',Color=grey(0.8) );
plot( mean(sf), mean(abs(zetaFloe).^2),'k-',LineWidth=1.5 );
ax.XScale = 'log';
ax.YScale = 'log';
ax.XLim = [50,50e3];
f = 1.4e-4;
ax.YLim = ([0.005,2]*f ).^2;

%% Save output

save('../MAT_Files/floeSpectra/idealSM_updated.mat','sf','UFloe','zetaFloe');
%}

%% ======== FUNCTIONS =====================================

function floeK = floeInventory(G,sf,floeResolution)
    if nargin<3, floeResolution = 150; end
%     minVerts = 8;
    minVerts = 24;
    M = length(G);
    m = randi(M,1);

    T = size(G{m},1);
    N = max( [ round(5*(sf/floeResolution)), minVerts+1] );
    t = 1:T; tp = linspace(1,T,N);
    fx = interp1(t,G{m}(:,1),tp).';
    fy = interp1(t,G{m}(:,2),tp).';
    pgon = polyshape( fx(1:N-1), fy(1:N-1), KeepCollinearPoints=true );

    [fcx,fcy] = centroid(pgon);
    floeArea = area(pgon);

    F(:,1) = fx(1:N-1)-fcx;
    F(:,2) = fy(1:N-1)-fcy;
    
    floeK = F * sqrt( (sf.^2)/floeArea );
end

