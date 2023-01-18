function floeOut = readFloeOut( outName, relativeOutPath, opts )
% READFLOEOUT reads data from FloeDyn output file
%
%   out = readFloeOut( outName )
%   out = readFloeOut( outName, relativeOutPath ) where 'relativeOutPath'
%   is the file path RELATIVE to '/Floe_Cpp/io/outputs/'
%
%   The output structure is organized as:
%       floeOut.fname           : FloeDyn filename
%       floeOut.fnameOriginal   : FloeDyn original filename      (optional)
%       floeOut.forceName       : Filename of forcing file used  (optional)
%       floeOut.processDate     : Date that this processing was run (yyyy-mm-dd)
%       floeOut.chunk:
%           floeOut.chunk.ke    : Kinetic energy
%           floeOut.chunk.obls  : Ocean Boundary Layer speed 
%           floeOut.chunk.mc    : Mass center
%           floeOut.chunk.win   : Domain window
%           floeOut.chunk.sic   : Domain sea ice concentration (initial)
%       floeOut.floe:
%           floeOut.floes.numFloes
%           floeOut.floes.floeSize
%           floeOut.floes.floeArea
%           floeOut.floes.floeShapes
%           floeOut.floes.floeMesh
%       floeOut.state:
%           floeOut.state.time : elapsed time [seconds]
%           floeOut.state.xa   : absolute x position [m]
%           floeOut.state.ya   : absolute y position [m]
%           floeOut.state.th   : rotation angle 
%           floeOut.state.ui   : zonal velocity [m/s]
%           floeOut.state.vi   : meridional velocity [m/s]
%           floeOut.state.om   : angular velocity 
%           floeOut.state.im   : total_received_impulse
%           floeOut.state.xp   : 'image' x position (for periodic BCs) [m]
%           floeOut.state.yp   : 'image' x position (for periodic BCs) [m]
%           floeOut.state.PBC  : booleen flag for periodic bcs
%
%   S.D.Brenner, 2022


%% Parse inputs/define directories

arguments
    outName             {mustBeText}
    relativeOutPath     {mustBeText} = '';
    opts.rootDir        {mustBeText} = '~/Documents/Brown/Projects/SASIP/FloeDyn/Floe_Cpp/io/outputs/';
    opts.orginalOutName {mustBeText} = '';
    opts.forceName      {mustBeText} = '';
end


%% Check for output file

fName = fullfile(opts.rootDir,relativeOutPath,outName);
if ~exist(fName,'file') 
    error('Could not find file %s in "%s"',...
        fullfile(relativeOutPath,outName),opts.rootDir );
end

                  
%% Create structure and add file/processing details 

floeOut.fname = fullfile(relativeOutPath,outName);
if ~isempty(opts.orginalOutName)
    floeOut.fnameOriginal = fullfile(relativeOutPath,opts.orginalOutName);
end
if ~isempty(opts.forceName)
    floeOut.forceName = opts.forceName;
end
floeOut.processDate = datestr( now(),'yyyy-mm-dd' );


%% Extract outputs from h5 file    


hfo = h5info(fName);
% Read in data from '/'
floeOut.chunk.ke   = h5read(fName,'/Kinetic Energy').';
floeOut.chunk.obls = h5read(fName,'/OBL_speed');
floeOut.chunk.mc   = h5read(fName,'/mass_center');
floeOut.state.time = h5read(fName,'/time'); % time in seconds
floeOut.chunk.win  = h5read(fName,'/window');

% Read in data from '/floe_shapes/' and '/floe_mesh/'
floeList = {hfo.Groups(1).Datasets.Name};
[~,sdx] = sort(str2double(floeList));
floeList = floeList( sdx );
M = length(floeList);
floeOut.floes.numFloes = M;
floeOut.floes.floeShapes = cell(1,M); 
floeOut.floes.floeMesh = cell(1,M); 
for m = 1:M
    floeOut.floes.floeShapes{m} = h5read(fName, ['/floe_shapes/',floeList{m}] ); 
    floeOut.floes.floeMesh{m} = h5read(fName, ['/floe_shapes/',floeList{m}] );  
end

%% Parse outputs

% 'Floe states'
fs = h5read(fName,'/floe_states');
[~,M,L] = size(fs);
floeOut.state.time = reshape( floeOut.state.time, 1,L );
floeOut.state.xa = reshape( fs(1,:,:), M,L );  % absolute x position
floeOut.state.ya = reshape( fs(2,:,:), M,L );  % absolute y position
floeOut.state.th = reshape( fs(3,:,:), M,L );  % rotation angle
floeOut.state.ui = reshape( fs(4,:,:), M,L );  % meridional velocity
floeOut.state.vi = reshape( fs(5,:,:), M,L );  % zonal velocity
floeOut.state.om = reshape( fs(6,:,:), M,L );  % angular velocity
floeOut.state.im = reshape( fs(7,:,:), M,L );  % total_received_impulse
floeOut.state.xp = reshape( fs(8,:,:), M,L );  % 'image' x position (for periodic bcs)
floeOut.state.yp = reshape( fs(9,:,:), M,L );  % 'image' x position (for periodic bcs)

% Check for periodic boundary conditions
if isequal(floeOut.state.xa,floeOut.state.xp) && isequal(floeOut.state.ya,floeOut.state.yp)
    floeOut.state.PBC = false;
else
    floeOut.state.PBC = true;
end

% Floe sizes
floeOut.floes.floeArea  = cellfun( @(C) polyarea(C(1,:),C(2,:)), floeOut.floes.floeShapes);
floeOut.floes.floeSize = sqrt(floeOut.floes.floeArea);

% Sea ice concentration
W = range( floeOut.chunk.win(1:2) );
H = range( floeOut.chunk.win(3:4) );
floeOut.chunk.sic = sum(floeOut.floes.floeArea) / ( W*H );

end