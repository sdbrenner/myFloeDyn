%% FLOE INITIALIZATION
% Generates a FloeDyn input file for a desired window/boundary area and sea
% ice distribution properties
%
% To use: adjust desired window/boundary properties in §1, floe-size
% distribution & sea ice targets in §2, and directory paths in §4 & §7,
% then RUN
%
% Notes: 
% - "floe size" is defined as the square root of the floe area
% - if a sea ice concentration (sic) target is specified, then the min/max
%   sizes of floes may be adjusted slightly to achieve the target perfectly
% - to achieve high sic (>75%), it is usually necessary to have a wide size
%   range (so that small floes can fit in between larger floes)
% - floe placement is by trial-and-error, and can be slow (especially for
%   high concentrations and high numbers of floes) 
% - for more than ~1000-2000 floes, plotting in §6 is pretty inefficient
% - written in MATLAB 2021b; relies on function 'getFloeSizes'.
%
% S.D.Brenner, 2022 

%% Clean workspace
clearvars;
clc;
close all;

%% §A1. Define window and boundary constraints
% The 'window' is included in the FloeDyn input file and defines the intial
% model domain. It is defined only be x and y extents, so is always a
% rectangle
%
% The 'boundary' can be any polygonal shape and defines the area in which
% the floes are distributed and the sea ice concentration is calculated.

% Define window
W = 50e3; H = 50e3;
X = 0; Y = 0; 
position = [X,Y,W,H]; %(for plotting)
window = [X+W*[0,1],Y+H*[0,1]];

% Define boundary:
    % boundary is the same as model domain window:
    boundary = [ X+[0.65;0.65;1;1;0.65]*W , Y+[0;1;1;0;0]*H];

    % boundary is a rectangle taking up half of the model domain:
    % boundary = [ X+[0.5;0.5;1;1;0.5]*W , Y+[0;1;1;0;0]*H];
    
    % boundary is a circle:
    % RT = W/2*exp(1i*linspace(-pi,pi,181))';
    % boundary = [real(RT),imag(RT)];

% Calculate area enclosed by boundary (used for calculating ice concentration)
Aw = polyarea( boundary(:,1), boundary(:,2) ); % boundary size

%% §A2. Define FSD and sea ice targets
% Ice concentration and floe number targets CANNOT be simulatenously
% specified for a prescribed FSD and 'boundary' area.
%
% Set either one of 'sicTarget' or 'numFloes', and leave the other as an
% empty vector. If values are given for both 'sicTarget' AND 'numFloes',
% then the script will shift the min and max sizes of the FSD in order to
% accomodate, which may be undesirable (and possibly unphysical)


% Define FSD properties
min_size = 1000;     % minimum floe size
max_size = 5000;   % maximum floe size
alpha = 2;          % non-cumulative FSD power-law exponent

% Targets: set either sea ice concentration ('sicTarget') or number of floes ('numFloes')
sicTarget = [0.5];
numFloes = [];


%% §A3. Generate and plot FSD


[floeSizes,ffun,Ffun] = getFloeSizes( [min_size,max_size,alpha],numFloes,sicTarget*Aw );
numFloes = length(floeSizes);
floeAreas = floeSizes.^2;
% floeAreas = sort(floeAreas,"descend");
sic = sum(floeAreas)/Aw;
fprintf('FSD generated with %g floes from %2.1f m to %2.1f m and ice concentration %2.1f%%\n',...
         numFloes, min(floeSizes), max(floeSizes), 100*sic )


fH = figure(1); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [8,3.5];
subplot(1,2,1);
hold on;

numBins = 20;
bins = logspace( log10(min_size),log10(max_size),numBins );
numFloesBin = histcounts(floeSizes,bins,'normalization','pdf'); 

loglog( bins, ffun(bins),'k-' );
loglog( binCenters(bins),numFloesBin,'o','color',lines(1),'markersize',4);
ax = gca; ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [ 0.9*min_size, 1.1*max_size ];
xlabel('floe size [m]');
ylabel('PDF');

subplot(1,2,2);
hold on;
% fs = sort(floeSizes,'ascend');
F = 1-(1:numFloes)/numFloes;
x = logspace( log10(min_size),log10(max_size),1e3 );
loglog( sort(floeSizes,'ascend'), F,'o','color',lines(1),'markersize',4);
loglog( x(1:end-1), Ffun(x(1:end-1)), 'k-','linewidth',1.5);
ax = gca; ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [ 0.9*min_size, 1.1*max_size ];
ax.YLim(2) = 5;
xlabel('floe size [m]');
ylabel('complementary CDF');

sgtitle('Floe Size Distribution')
drawnow;

%% §A4. Floe shape inventory

% Load default FloeDyn inventory of floe shapes 
load('../Floe_Cpp/io/inputs/Biblio_Floes.mat');

% Use idealized circular floes
% G{1} = getCircularFloe();



%% §A5. Generate and place floes
% Floes are randomly selected from the inventory, then randomly rotated and
% placed within the boundary starting from the largest and working down the
% list. After each placement, the floe position is rejected if it overlaps
% with another floe or the edge of the boundary. To avoid infinite loops, a
% limit is placed on the number of "rejections" allowed for any floe--if
% the limit is reached then placement fails and a new attempt is made from
% scratch. Convergence is not guaranteed, and is generally difficult for a
% high concentrations.


% Initialize constraints and counts
rejectCount = zeros(1,numFloes);
rejectThreshold = 5e3;
attemptCount = 1;
attemptCountThreshold = 3;
floe_states = zeros(9,numFloes);

% Pre-allocate other variables
clear floePos;
floePos = cell(1,numFloes);
floeList = cell(1,numFloes);
floeXY = cell(1,numFloes);

% Loop through remaning floes
disp('placing floes...')
k = 1;
while k<=numFloes

    % Get random floe from inventory
    floeK = floeInventory( G,floeSizes(k),50 );
    % Generate random floe position and orientation
    theta = 2*pi*rand(1);
    fsXY(1) = W*rand(1)+X;
    fsXY(2) = H*rand(1)+Y;
    R = [ cos(theta), -sin(theta) ; 
          sin(theta),  cos(theta)  ];
    % Put floe into position
    floePos{k} = (floeK*R)+fsXY;

    % Create variables used for making the HDF5 file
    floeList{k} = sprintf('%u',k-1);
    floeXY{k} = (floeK*R)';
    floe_states(1:2,k) = fsXY;
    floe_states(8:9,k) = fsXY;

    % Check that the floe doesn't overlap with any other floes
    [hitTestResult] = hitTest( floePos{k},floePos(1:k-1),boundary );
    if ~hitTestResult % accept and advance loop
        k = k+1;
    else % try again
        if rejectCount(k)>rejectThreshold
            warning('on','all');
            warning('Could not place all floes on attempt %g',attemptCount);
            attemptCount = attemptCount+1;
            if attemptCount>=attemptCountThreshold
                break;
            end 
            rejectCount = zeros(1,numFloes);
            k = 1;  
        end
        rejectCount(k)= rejectCount(k)+1;
    end

end
numPlacements = k-1;

% 
if numPlacements==numFloes
    disp('floe placement succeeded');
elseif attemptCount>=attemptCountThreshold
    warning('Could not place all floes after %g attempts',attemptCount);
end


% save inital FSD results
min_sizeA = min_size;
max_sizeA = max_size;
floeSizesA = floeSizes; 
numFloesA = numFloes;

%% §B1. Define secondary window and boundary constraints
% The 'window' is included in the FloeDyn input file and defines the intial
% model domain. It is defined only be x and y extents, so is always a
% rectangle
%
% The 'boundary' can be any polygonal shape and defines the area in which
% the floes are distributed and the sea ice concentration is calculated.

% Define boundary:
    % boundary is the same as model domain window:
    boundary = [ X+[0.55;0.55;0.75;0.75;0.55]*W , Y+[0;1;1;0;0]*H];

    % boundary is a rectangle taking up half of the model domain:
    % boundary = [ X+[0.5;0.5;1;1;0.5]*W , Y+[0;1;1;0;0]*H];
    
    % boundary is a circle:
    % RT = W/2*exp(1i*linspace(-pi,pi,181))';
    % boundary = [real(RT),imag(RT)];

% Calculate area enclosed by boundary (used for calculating ice concentration)
Aw = polyarea( boundary(:,1), boundary(:,2) ); % boundary size

%% §B2. Define secondary FSD and sea ice targets
% Ice concentration and floe number targets CANNOT be simulatenously
% specified for a prescribed FSD and 'boundary' area.
%
% Set either one of 'sicTarget' or 'numFloes', and leave the other as an
% empty vector. If values are given for both 'sicTarget' AND 'numFloes',
% then the script will shift the min and max sizes of the FSD in order to
% accomodate, which may be undesirable (and possibly unphysical)


% Define FSD properties
min_size = 150;     % minimum floe size
max_size = 1500;   % maximum floe size
alpha = 2;          % non-cumulative FSD power-law exponent

% Targets: set either sea ice concentration ('sicTarget') or number of floes ('numFloes')
sicTarget = [0.25];
numFloes = [];


%% §B3. Generate and plot secondary FSD

[floeSizes,ffun,Ffun] = getFloeSizes( [min_size,max_size,alpha],numFloes,sicTarget*Aw );
numFloes = length(floeSizes);
floeAreas = floeSizes.^2;
% floeAreas = sort(floeAreas,"descend");
sic = sum(floeAreas)/Aw;
fprintf('FSD generated with %g floes from %2.1f m to %2.1f m and ice concentration %2.1f%%\n',...
         numFloes, min(floeSizes), max(floeSizes), 100*sic )


numFloesAll = numFloes + numFloesA;
floeSizesAll = [floeSizesA,floeSizes];

fH = figure(1); 
fH.Units = 'inches';
fH.Position([3,4]) = [8,3.5];

subplot(1,2,1);
hold on;

numBins = 20;
bins = logspace( log10(min_size),log10(max_sizeA),numBins );
numFloesBin = histcounts(floeSizesAll,bins,'normalization','pdf'); 

loglog( bins, ffun(bins),'k-' );
loglog( binCenters(bins),numFloesBin,'o','markersize',4);
ax = gca; ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [ 0.9*min_size, 1.1*max_sizeA ];
xlabel('floe size [m]');
ylabel('PDF');

subplot(1,2,2);
hold on;
% fs = sort(floeSizes,'ascend');
F = 1-(1:numFloesAll)/numFloesAll;
x = logspace( log10(min_size),log10(max_size),1e3 );
loglog( sort(floeSizesAll,'ascend'), F,'o','markersize',4);
loglog( x(1:end-1), Ffun(x(1:end-1)), 'k-','linewidth',1.5);
ax = gca; ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [ 0.9*min_size, 1.1*max_sizeA ];
ax.YLim(2) = 5;
xlabel('floe size [m]');
ylabel('complementary CDF');

sgtitle('Floe Size Distribution')
drawnow;



%% §B5. Generate and place floes
% Floes are randomly selected from the inventory, then randomly rotated and
% placed within the boundary starting from the largest and working down the
% list. After each placement, the floe position is rejected if it overlaps
% with another floe or the edge of the boundary. To avoid infinite loops, a
% limit is placed on the number of "rejections" allowed for any floe--if
% the limit is reached then placement fails and a new attempt is made from
% scratch. Convergence is not guaranteed, and is generally difficult for a
% high concentrations.


% Initialize constraints and counts
rejectCount = zeros(1,numFloes);
rejectThreshold = 5e3;
attemptCount = 1;
attemptCountThreshold = 3;


% Loop through remaning floes
disp('placing floes...')
k = numFloesA+1;
kB = 1;
while kB<=numFloes

    % Get random floe from inventory
    floeK = floeInventory( G,floeSizes(kB),50 );
    % Generate random floe position and orientation
    theta = 2*pi*rand(1);
    fsXY(1) = W*rand(1)+X;
    fsXY(2) = H*rand(1)+Y;
    R = [ cos(theta), -sin(theta) ; 
          sin(theta),  cos(theta)  ];
    % Put floe into position
    floePos{k} = (floeK*R)+fsXY;

    % Create variables used for making the HDF5 file
    floeList{k} = sprintf('%u',k-1);
    floeXY{k} = (floeK*R)';
    floe_states(1:2,k) = fsXY;
    floe_states(8:9,k) = fsXY;

    % Check that the floe doesn't overlap with any other floes
    [hitTestResult] = hitTest( floePos{k},floePos(1:k-1),boundary );
    if ~hitTestResult % accept and advance loop
        k = k+1;
        kB = kB+1;
    else % try again
        if rejectCount(kB)>rejectThreshold
            warning('on','all');
            warning('Could not place all floes on attempt %g',attemptCount);
            attemptCount = attemptCount+1;
            if attemptCount>=attemptCountThreshold
                break;
            end 
            rejectCount = zeros(1,numFloes);
            k = numFloesA+1;
            kB = 1;
        end
        rejectCount(kB)= rejectCount(kB)+1;
    end

end
numPlacements = k-1;

% 
if numPlacements==numFloes
    disp('floe placement succeeded');
elseif attemptCount>=attemptCountThreshold
    warning('Could not place all floes after %g attempts',attemptCount);
end

%% §C1. Plot


%{
% Plot trials    
rejectCount(1:numFloes>numPlacements+1) = [];

fH = figure(2); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [8,3.5];

subplot(1,2,1); 
semilogy( rejectCount,'ko'); 
yline(rejectThreshold,'k',LineWidth=1);
xline(numPlacements,'r--',LineWidth=1);
grid on; 
xlim([1,numFloes])
xlabel('floe #');
ylabel('placement attempts');

subplot(1,2,2); 
semilogy( cumsum(rejectCount),'ko'); 
xline(numPlacements,'r--',LineWidth=1);
grid on; 
xlim([1,numFloes])
xlabel('floe #');
ylabel('cumulative placement attempts');

sgtitle('Floe Placement Attempts')

%}

% Plot floe positions
fH = figure(10); clf; hold on;

aspect = position(3)/(position(4));
fH.Units = 'inches';
fH.Position(3) = aspect*fH.Position(4);

rH = rectangle('Position',position);
rH.LineStyle = ':';

% plot( boundary(:,1), boundary(:,2),'r--');


ax = gca;
ax.XLim = X + [0,1]*W+0.1*[-1,1]*W;
ax.YLim = Y + [0,1]*H+0.1*[-1,1]*W;
daspect([1,1,1]);
drawnow;
for k = 1:numPlacements
    figure(10);
    patch( floePos{k}(:,1), floePos{k}(:,2), grey(0.9) );
%     drawnow;
end
% patchfun = @(C) patch( C(:,1), C(:,2), grey(0.95),'edgecolor',grey(0.35) );
% cellfun( patchfun, floePos );
% drawnow; pause(0.5);



if attemptCount>=attemptCountThreshold
%     error('Could not place all floes after %g attempts',attemptCount);
    return;
end




%% §C2. Create HDF5 file

% Generate automatic filename
if log10(W)>4 && log10(H)>4
    inName = sprintf("in_%2.0fp_%ix%ikm_%if",100*sic,round(W/1e3),round(H/1e3),numFloes);
else
    inName = sprintf("in_%2.0fp_%ix%im_%if",100*sic,round(W),round(H),numFloes);
end
inName = strrep(inName,' ','');  % remove spaces in name
inName = strrep(inName,'.',','); % remove periods in name



prompt = sprintf('save hdf5 input file (''%s''),   y/n/filename [y]?: ',inName);
inputStr = input(prompt,'s');
if isempty(inputStr)
    inputStr = 'y';
end
% if 'n' or 'no', exit script (don't save) 
if any( strcmpi( inputStr,{'n','no'}) ) 
    return;
% if NOT 'n' or 'no', AND NOT 'y' or 'yes', then a new file-name is given.
% Overwrite inName:
elseif ~strcmpi( inputStr, {'y','yes'})
    inName = inputStr;
end
% if the file extension isn't already given, add it:
inName = convertStringsToChars(inName);
if ~strcmp(inName(end-2:end),'.h5')
    inName = strcat(inName,".h5");   % add file extension 
end

fPath = "../Floe_Cpp/io/inputs/";
%     fPath = "../Floe_Cpp/io/inputs/50x50km_N_specified/";
fName = strcat(fPath,inName);
if exist(fName,'file'), delete(fName); end

% Create empty HDF5 file
h5create(fName,'/floe_states',size(floe_states,1:3),'ChunkSize',size(floe_states,1:3));
h5create(fName,'/window',length(window) );
for k = 1:length(floeList)
%     S = hfo.Groups.Datasets(k).Dataspace.Size;
    S = size(floeXY{k});
    h5create(fName, ['/floe_shapes/',floeList{k}], S );
end
% Write data to file
h5write(fName,'/floe_states',floe_states);
h5write(fName,'/window',window);
for k = 1:length(floeList)
    % check overwrite size
%     S = hfo.Groups.Datasets(k).Dataspace.Size(2);
    h5write(fName, ['/floe_shapes/',floeList{k}], floeXY{k} );
end
fprintf('Saved File: %s\n',fName)




%% §D. Embedded functions


function floeK = floeInventory(G,sf,floeResolution)
    if nargin<3, floeResolution = 150; end
    minVerts = 8; % minimum number of vertices to describe floe
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




function [hitTestResult] = hitTest( floePos,allFloes,boundary )
    % chech if any points are outside window
    hitTestResult = any( ~inpolygon( ...
            floePos(:,1),floePos(:,2), ...
            boundary(:,1),boundary(:,2)...
            ));
    % check if any points are inside other floes
    if ~isempty(allFloes)
        for k = 1:length(allFloes)
            p1x = floePos(:,1);
            p1y = floePos(:,2);
            p2x = allFloes{k}(:,1);
            p2y = allFloes{k}(:,2);
            if any( [inpolygon(p1x,p1y,p2x,p2y) ; inpolygon(p2x,p2y,p1x,p1y)] )
                hitTestResult = 1;
                break;
            end
        end
    end

end



function  floeXY = getCircularFloe()

    theta = linspace(0,2*pi);
    [x,y] = pol2cart(theta,100);
    floeXY = [x(:),y(:)];
end



