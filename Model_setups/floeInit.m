%% FLOE INITIALIZATION
% Generates a FloeDyn input file for a desired window/boundary area and sea
% ice distribution properties
%
% To use, adjust:
%   (a) desired window/boundary properties in §1; 
%   (b) floe-size distribution & sea ice targets in §2; and,
%   (c) directory paths in §4 & §7;
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
% - Redraws (interpolates) floe boundaries with a set resolution and
%   minimum number of vertices (to adjust, see 'floeInventory' function).    
% - for more than ~1000-2000 floes, plotting in §6 is pretty inefficient
% - written in MATLAB 2021b; relies on function 'getFloeSizes'.
%
% S.D.Brenner, 2022 

%% Clean workspace
clearvars;
clc;
% close all;

%% §1. Define window and boundary constraints
% The 'window' is included in the FloeDyn input file and defines the intial
% model domain. It is defined only by x and y extents, so is always a
% rectangle
%
% The 'boundary' can be any polygonal shape and defines the area in which
% the floes are distributed and the sea ice concentration is calculated.

% Define window
% W = 500e3; H = 500e3;
W = 160e3; H = 90e3;
X = 0; Y = 0; 
position = [X,Y,W,H]; %(for plotting)
window = [X+W*[0,1],Y+H*[0,1]];

% Define boundary:
    % boundary is the same as model domain window:
    boundary = [ X+[0;0;1;1;0]*W , Y+[0;1;1;0;0]*H];

    % boundary is a rectangle taking up half of the model domain:
%     boundary = [ X+[0;0;1;1;0]*W , Y+[0.5;1;1;0.5;0.5]*H];
    
    % boundary is a circle:
    % RT = W/2*exp(1i*linspace(-pi,pi,181))';
    % boundary = [real(RT),imag(RT)];

% Calculate area enclosed by boundary (used for calculating ice concentration)
Abound = polyarea( boundary(:,1), boundary(:,2) ); % boundary size
Awin = W*H; % window size

%% §2. Define FSD and sea ice targets
% Ice concentration and floe number targets CANNOT be simulatenously
% specified for a prescribed FSD and 'boundary' area.
%
% Set either one of 'sicTarget' or 'numFloes', and leave the other as an
% empty vector. If values are given for both 'sicTarget' AND 'numFloes',
% then the script will shift the min and max sizes of the FSD in order to
% accomodate, which may be undesirable 


% Define FSD properties
% minSize = 3e3;     % minimum floe size
% maxSize = 30e3;   % maximum floe size
minSize = 300;     % minimum floe size
maxSize = 5000;   % maximum floe size
alpha = 2;          % non-cumulative FSD power-law exponent
sizeSpan = (maxSize-minSize)./(maxSize+minSize);


% Targets: set either sea ice concentration ('sicTarget') or number of floes ('numFloes')
sicTarget = 0.65;
numFloes = [];


%% §3. Generate and plot FSD


floeSizes = getFloeSizes( [minSize,maxSize,alpha],numFloes,sicTarget*Abound );
% minSize = min(floeSizes); maxSize = max(floeSizes);
[ffun,Ffun] = fsdFuns( minSize,maxSize,alpha ); 

numFloes = length(floeSizes);
floeAreas = floeSizes.^2;
% floeAreas = sort(floeAreas,"descend");
sicBound = sum(floeAreas)/Abound;
sicWin = sum(floeAreas)/Awin;
% Displace message with FSD/SIC properties
if abs(sicWin-sicBound)<0.01 % if boundary and windown are the same sic
    fprintf(['FSD generated with %g floes from %2.1f m to %2.1f m\n',...
             'sea ice concentration %2.1f%%\n\n'],...
             numFloes, min(floeSizes), max(floeSizes), 100*sicBound );
else
    fprintf(['FSD generated with %g floes from %2.1f m to %2.1f m\n',...
                 'sea ice concentration %2.1f%% within the boundary\n'...
                 '(%2.1f%% within the window)\n\n'],...
                 numFloes, min(floeSizes), max(floeSizes), 100*sicBound,100*sicWin );    
end

fH = figure(1); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [8,3.5];
subplot(1,2,1);
hold on;

% Design bins
numBins = 20;
% Logarithmic spacing:
bins = logspace( log10(minSize),log10(maxSize),numBins );
% Spacing based on CCDF:
% x = ( logspace( log10(min_size),log10(max_size),1e3 ) );
% bins = interp1(Ffun(x),x,linspace(0,1,numBins));
% bins(end) = x(1);
% bins = sort(bins);

numFloesBin = histcounts(floeSizes,bins,'normalization','pdf'); 

loglog( bins, ffun(bins),'k-' );
loglog( binCenters(bins),numFloesBin,'o','color',lines(1),'markersize',4);
ax = gca; ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [ 0.9*minSize, 1.1*maxSize ];
xlabel('floe size [m]');
ylabel('PDF');
YL = ylim;
plot( bins, YL(2),'|',Color=0.7*[1,1,1]); 
ax.YLim = YL;

subplot(1,2,2);
hold on;
% fs = sort(floeSizes,'ascend');
F = 1-(1:numFloes)/numFloes;
x = logspace( log10(minSize),log10(maxSize),1e3 );
PL = 1.2*(x./minSize).^(1-alpha);
loglog( sort(floeSizes,'ascend'), F,'o','color',lines(1),'markersize',4);
loglog( x(1:end-1), Ffun(x(1:end-1)), 'k-','linewidth',1.5);
loglog( x, PL, '--',Color=[0.5,0.5,0.5] )
ax = gca; ax.XScale = 'log'; ax.YScale = 'log';
ax.XLim = [ 0.9*minSize, 1.1*maxSize ];
ax.YLim(2) = 5;
xTxt = minSize + 0.1*maxSize;
yTxt = 1.6*(xTxt./minSize).^(1-alpha);
text( xTxt, yTxt,sprintf('$\\sim{s_f^{%2.2g}}$',1-alpha),Color=[0.5,0.5,0.5],Interpreter="latex" );
xlabel('floe size [m]');
ylabel('complementary CDF');

sgtitle('Floe Size Distribution')
drawnow;

%% §4. Floe shape inventory

% Load default FloeDyn inventory of floe shapes 
load('../../Floe_Cpp/io/inputs/Biblio_Floes.mat');

% Use idealized circular floes
% G{1} = getCircularFloe();



%% §5. Generate and place floes
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
updateLabIncrement = 0.1;
updateLab = updateLabIncrement;
while k<=numFloes

    % Get random floe from inventory
%     floeResolution = 150; % for submesoscale sim
%     floeResolution = 1500;  % for 2DQG sim
    floeResolution = minSize/2;
    floeK = floeInventory( G,floeSizes(k),floeResolution );

    % Generate random floe position and orientation
    theta = 2*pi*rand(1);
    fsXY(1) = W*rand(1)+X;
    fsXY(2) = H*rand(1)+Y;
    R = [ cos(theta), -sin(theta) ; 
          sin(theta),  cos(theta)  ];
    % Put floe into position
    floePos{k} = (floeK*R)+fsXY;

    % Create variables used for making the HDF5 file
%     floeList{k} = sprintf('%u',k-1);
    floeXY{k} = (floeK*R)';
    floe_states(1:2,k) = fsXY;
    floe_states(8:9,k) = fsXY;

    % Check that the floe doesn't overlap with any other floes
    [hitTestResult] = hitTest( floePos{k},floePos(1:k-1),boundary );
    if ~hitTestResult % accept and advance loop
        k = k+1;
        if (k/numFloes)>=updateLab
            fprintf('%02g%% placed...\n',100*updateLab );
            updateLab = updateLab+updateLabIncrement;
        end
    else % try again
        if rejectCount(k)>rejectThreshold
%             warning('on','all');
            warning('Could not place all floes on attempt %g (%02g%% placed)',...
                    attemptCount, 100*(k/numFloes) );
            % reset initialization
            rejectCount = zeros(1,numFloes);
            k = 1;  
            updateLab = updateLabIncrement;
            floePos = cell(1,numFloes);
            floeXY = cell(1,numFloes);
            % update attempt counter
            attemptCount = attemptCount+1;
            if attemptCount>attemptCountThreshold
                break;
            end 

        end
        rejectCount(k)= rejectCount(k)+1;
    end

end
numPlacements = k-1;

 
if numPlacements==numFloes
    fprintf('floe placement succeeded\n\n');
elseif attemptCount>attemptCountThreshold
    warning('Could not place all floes after %g attempts\n',attemptCountThreshold);
end

%% §6. Plot


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
figNum = 100;
fH = figure(figNum); clf; hold on;

aspect = position(3)/(position(4));
fH.Units = 'inches';
fH.Position(3) = aspect*fH.Position(4);

rH = rectangle('Position',position);
rH.LineStyle = ':';

% plot( boundary(:,1), boundary(:,2),'y--');


ax = gca;
ax.XLim = X + [0,1]*W+0.05*[-1,1]*W;
ax.YLim = Y + [0,1]*H+0.05*[-1,1]*W;
ax.Color = [30,150,220]/255;
daspect([1,1,1]);
% for k = 1:numPlacements
%     figure(figNum);
%     patch( floePos{k}(:,1), floePos{k}(:,2), [0.95,0.95,0.9],...
%            'edgecolor','none');
% %     drawnow;
% end
patchfun = @(C) patch( C(:,1), C(:,2), [0.95,0.95,0.98],'edgecolor','none' );
floePosPlot = floePos( cellfun( @(C) numel(C)>0, floePos ) );
cellfun( patchfun, floePosPlot );


%% §7. Create HDF5 file

% Don't try to save if floes weren't all placed
if attemptCount>attemptCountThreshold
    return;
end


% GENERATE AUTOMATIC FILENAME
if log10(W)>4 && log10(H)>4
    inName = sprintf("in_%02.0fp_%ix%ikm_%if",100*sicBound,round(W/1e3),round(H/1e3),numFloes);
else
    inName = sprintf("in_%02.0fp_%ix%im_%",if100*sicBound,round(W),round(H),numFloes);
end
inName = strrep(inName,' ','');  % remove spaces in name
inName = strrep(inName,'.',','); % remove periods in name (replace with commas)

% SAVE PROMPT
prompt = sprintf('save hdf5 input file (''%s'')?\n   y/n/filename [y]: ',inName);
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


% FILE PATH
mainFPath = "../../Floe_Cpp/io/inputs/";
% fSubDir = sprintf("sic_batches/FSD_%g_%g_%g/",minSize,maxSize,alpha);
% fSubDir = "sic_batches/FSD_300_7500_2/";
fSubDir = "custom/";
fPath = strcat( mainFPath,fSubDir );
% Check that path exists
if ~exist(fPath,"dir")
    warning('File path %s does not exist',fPath);
    prompt = sprintf('Make directory?\n   y/n [y]: ');
    inputStr = input(prompt,'s');
    if isempty(inputStr) || any( strcmpi( inputStr,{'y','yes'}) )
        mkdir(fPath);
    end
end


% SAVE FILE
fName = strcat(fPath,inName);
if exist(fName,'file'), delete(fName); end
% Create empty HDF5 file
h5create(fName,'/floe_states',size(floe_states,1:3),'ChunkSize',size(floe_states,1:3));
h5create(fName,'/window',length(window) );
for k = 1:length(floeList)
%     S = hfo.Groups.Datasets(k).Dataspace.Size;
    S = size(floeXY{k});
    floeLabel = sprintf('%u',k-1);
    h5create(fName, ['/floe_shapes/',floeLabel], S );
end
% Write data to file
h5write(fName,'/floe_states',floe_states);
h5write(fName,'/window',window);
for k = 1:length(floeXY)
    % check overwrite size
%     S = hfo.Groups.Datasets(k).Dataspace.Size(2);
    floeLabel = sprintf('%u',k-1);
    h5write(fName, ['/floe_shapes/',floeLabel], floeXY{k} );
end
fprintf('Saved File: %s\n',fName)




%% §8. Embedded functions



function floeK = floeInventory(G,sf,floeResolution)
    if nargin<3, floeResolution = 150; end
    minVerts = 6; % minimum number of vertices to describe floe
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


% % overload input to ignore prompt:
% function out = input(varargin)
%     out = 'y';
% end
