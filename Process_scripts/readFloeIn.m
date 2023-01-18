function floeIn = readFloeIn( inName, relativePath, opts )
% READFLOEIN reads data from FloeDyn input file
%
%   out = readFloeIn( outName )
%   out = readFloeIn( outName, relativeOutPath ) where 'relativeOutPath'
%   is the file path RELATIVE to '/Floe_Cpp/io/outputs/'
%
%   The output structure is organized as:
%       floeIn.fname           : FloeDyn filename
%       floeIn.chunk:
%           floeIn.chunk.win   : Domain window
%           floeIn.chunk.sic   : Domain sea ice concentration (initial)
%       floeIn.floe:
%           floeIn.floes.numFloes
%           floeIn.floes.floeSize
%           floeIn.floes.floeArea
%           floeIn.floes.floeShapes
%       floeIn.state:
%           floeIn.state.xa   : absolute x position [m]
%           floeIn.state.ya   : absolute y position [m]
%           floeIn.state.th   : rotation angle 
%           floeIn.state.ui   : zonal velocity [m/s]
%           floeIn.state.vi   : meridional velocity [m/s]
%           floeIn.state.om   : angular velocity 
%           floeIn.state.im   : total_received_impulse
%           floeIn.state.xp   : 'image' x position (for periodic BCs) [m]
%           floeIn.state.yp   : 'image' x position (for periodic BCs) [m]
%
%   S.D.Brenner, 2022


%% Parse inputs/define directories

arguments
    inName             {mustBeText}
    relativePath        {mustBeText} = '';
    opts.rootDir        {mustBeText} = '~/Documents/Brown/Projects/SASIP/FloeDyn/Floe_Cpp/io/inputs/';
end


%% Check for output file

fName = strcat(opts.rootDir,relativePath,inName);
if ~exist(fName,'file') 
    error('Could not find file %s in "%s"',...
        strcat(relativePath,inName),opts.rootDir );
end

                  
%% Create structure and add file/processing details 

floeIn.fname = strcat(relativePath,inName);



%% Extract outputs from h5 file    


hfo = h5info(fName);
% Read in data from '/'
floeIn.chunk.win  = h5read(fName,'/window');

% Read in data from '/floe_shapes/' and '/floe_mesh/'
floeList = {hfo.Groups(1).Datasets.Name};
[~,sdx] = sort(str2double(floeList));
floeList = floeList( sdx );
M = length(floeList);
floeIn.floes.numFloes = M;
floeIn.floes.floeShapes = cell(1,M); 
floeIn.floes.floeMesh = cell(1,M); 
for m = 1:M
    floeIn.floes.floeShapes{m} = h5read(fName, ['/floe_shapes/',floeList{m}] ); 
end

%% Parse inputs

% 'Floe states'
fs = h5read(fName,'/floe_states');
[~,M,L] = size(fs);
floeIn.state.xa = reshape( fs(1,:,:), M,L );  % absolute x position
floeIn.state.ya = reshape( fs(2,:,:), M,L );  % absolute y position
floeIn.state.th = reshape( fs(3,:,:), M,L );  % rotation angle
floeIn.state.ui = reshape( fs(4,:,:), M,L );  % meridional velocity
floeIn.state.vi = reshape( fs(5,:,:), M,L );  % zonal velocity
floeIn.state.om = reshape( fs(6,:,:), M,L );  % angular velocity
floeIn.state.im = reshape( fs(7,:,:), M,L );  % total_received_impulse
floeIn.state.xp = reshape( fs(8,:,:), M,L );  % 'image' x position (for periodic bcs)
floeIn.state.yp = reshape( fs(9,:,:), M,L );  % 'image' x position (for periodic bcs)


% Floe sizes
floeIn.floes.floeArea  = cellfun( @(C) polyarea(C(1,:),C(2,:)), floeIn.floes.floeShapes);
floeIn.floes.floeSize = sqrt(floeIn.floes.floeArea);

% Sea ice concentration
W = range( floeIn.chunk.win(1:2) );
H = range( floeIn.chunk.win(3:4) );
floeIn.chunk.sic = sum(floeIn.floes.floeArea) / ( W*H );

end