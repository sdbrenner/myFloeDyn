%% Animate floeDyn output file

%% Clean workspace
clear;
clc;
close all;

%% Update paths

addpath('../Process_scripts/');
addpath('../Analysis/');

%% Define floe and forcing files

% FloeDyn output
outPath = '';
% outPath = 'long_runs/';
% outName = 'out_30p_50x50km_1237f_geo.h5'; % LONG
% outName = 'out_QVdTK.h5'; % Collision!
% outName = 'out_78wt9.h5';
% outName = 'out_IIIPQ.h5';
% outName = 'out_AKbB0.h5'; % long run no geo
% outName = 'out_0UMab.h5';
% outName = 'out_uG2LR.h5';
% outName = 'out_EoopZ.h5';o

outName = 'out_zOBFG.h5';


% forcing
forcePath = 'forcings/';
forceName = [];
% forceName = '_input_forcing_SMgeo_long.mat';
% forceName = '_input_forcing_SMgeo_5day.mat';
% forceName = '_input_forcing_2DQGn.mat';

%% Read in data


fo = readFloeOut( outName, outPath );
fo  = makeDomainPoly(fo);
L = length(fo.chunk.icePoly);

% if ~isempty(forceName)
%     forcing = readForcing(forceName,forcePath);
%     tol = 120; % timestep tolerance [seconds]
%     [LIA,~] = ismembertol( fo.state.time, forcing.t,tol,'DataScale',1 );
%     fo = subsetFloeOut( fo, LIA );
% else
    forcing = [];
% end

%% Read in forcing and interpolate to hourly values

% if ~exist('forcing','var')
%     forcing = readForcing(forceName,forcePath);
% end
% 
% 
% [fx,fy,ot] = meshgrid(forcing.x, forcing.y, fo.state.time);
% f2.x = forcing.x; 
% f2.y = forcing.y;
% f2.t = fo.state.time;
% 
% f2.u = interp3( forcing.x, forcing.y, forcing.t, forcing.u,...
%                 fx,fy,ot );
% f2.v = interp3( forcing.x, forcing.y, forcing.t, forcing.v,...
%                 fx,fy,ot );
% f2.u(:,:,end) = forcing.u(:,:,end);
% f2.v(:,:,end) = forcing.v(:,:,end);
%% Set figure/plotting options

% Figure/plotting options (general)
figOpts = struct;
vidFilePrefix = 'vid';
figOpts.aspectRatio = 2;
% figOpts.aspectRatio = 16/9;

% figOpts.iceBorder = grey(0.2);

% Background pcolor based on vorticity
% f = 1.4e-4;
% figOpts.iceColour = grey(0.95); 
% figOpts.iceBorder = grey(0.2);
% figOpts.iceAlpha  = 0.7;
% figOpts.cmap = cbrewer2('RdBu',21);
% figOpts.pFun = @(X,Y,U,V) 2*curl(X,Y,U,V)/f;  
% figOpts.CLim = 2*[-1,1];
% figOpts.cbLab = 'ocean surface vorticity (\zeta/f)';
% figOpts.axColour = 'none';
% vidFilePrefix = 'vrt';

% Background pcolor based on divergence
% figOpts.iceColour = 'none'; 
% figOpts.iceBorder = grey(0.2);
% figOpts.iceAlpha  = 0.7;
% figOpts.cmap = cbrewer2('RdBu',21); figOpts.cmap(11,:)=grey(1);
% figOpts.pFun = @(X,Y,U,V) divergence(X,Y,U,V); 
% figOpts.CLim = 1.55e-4*[-1,1];   
% figOpts.cbLab = 'ocean surface divergence [1/s]';
% figOpts.axColour = 'none';
% vidFilePrefix = 'div';


%% Make and animate figure

% Convert figure options to a cell of Name-Value pairs
figOptsNames = string(fieldnames(figOpts));
figOptsValues = struct2cell(figOpts);
figOptsNameValuePairs = {figOptsNames{:};figOptsValues{:}};

% Run animation
imFrame = animateFloes( fo,forcing,figOptsNameValuePairs{:} );
% imFrame = animateFloes( fo,f2,figOptsNameValuePairs{:} );
figPos = get(gcf,'Position');

%% Playback movie

imFrameS = cellstruct2structarray(imFrame);

fH = figure(10); 
fH.Position = figPos;
ax = gca;
ax.Position= [0,0,1,1];
movie(imFrameS,1)

%% Save video or gif

saveVideo = 1;
saveGif = 0;

if saveVideo
    vidPath = '../../Floe_Cpp/io/videos/vids/';
    vidName = [vidFilePrefix,outName(4:end-3),'_wide'];
    makeVid( [vidPath,vidName],imFrame,60 );
end


if saveGif
    gifPath = '../../Floe_Cpp/io/videos/gifs/';
    gifName = ['gif',outName(4:end-3)];
    makeGif( [gifPath,gifName,'.gif'], imFrame );
end


%%

function makeGif(fName,imFrame)
    
    delaytime = 1/8;
    % create gif file
    k = 1;
    K = length(imFrame);
    [A,map] = rgb2ind( imFrame{k}.cdata,256 );
    imwrite( A,map,fName,'gif',...
             LoopCount=Inf,DelayTime=delaytime);
    % append frames
    for k = 2:K
        [A,map] = rgb2ind( imFrame{k}.cdata,256 );
        imwrite( A,map,fName,'gif',...
                 WriteMode='append',DelayTime=delaytime);
    end
    disp('gif saved');
    fprintf('saved %s\n',fName);
end



function makeVid(fName,imFrame,videoLength)

    if nargin<3; videoLength = 30; end %[sec]
    pauseTime = 2; % freeze frame at the end of the video
    K = length(imFrame);
    fps = K/(videoLength-pauseTime);
    % Add a 3-second freeze frame at the end of the video (repeat last frame);
    for k = 1:ceil(pauseTime*fps)
        imFrame{K+k} = imFrame{K};
    end

    writerObj = VideoWriter(fName,'MPEG-4');
    writerObj.FrameRate = fps;
    writerObj.Quality = 100;

    open(writerObj);
    for k = 1:numel(imFrame)
        writeVideo(writerObj,imFrame{k} );
    end
    close(writerObj);
    fprintf('saved %s.%s\n',fName,writerObj.FileFormat);

   
end


function Y = cellstruct2structarray(X)
    N = length(X);
    Y(N) = X{N};
    for n = 1:length(X)
        Y(n) = X{n};
    end
end
