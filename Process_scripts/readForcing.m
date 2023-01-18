function  forcing = readForcing( forceName, relativeForcePath,opts )
% READFORCING reads data from FloeDyn input forcing file
%
%   forcing = readForcing( forceName );
%   forcing = readForcing( forceName, relativeForcePath ) where
%   'relativeForcePath' is the file path RELATIVE to '/Floe_Cpp/io/inputs/'
%
%   S.D.Brenner, 2022


%% Parse inputs/define directories

arguments
    forceName           {mustBeText}
    relativeForcePath   {mustBeText} = '';
    opts.rootDir        {mustBeText} = '~/Documents/Brown/Projects/SASIP/FloeDyn/Floe_Cpp/io/inputs/';
    opts.calcFields     {mustBeNumericOrLogical} = true;
end

%% Check for forcing file

fName = strcat(opts.rootDir,relativeForcePath,forceName);
if ~exist(fName,'file') 
    error('Could not find file %s in "Floe_Cpp/io/inputs/"',...
        strcat(relativeForcePath,forceName) );
end

%% Load forcing data into structure

forcing = load(fName);

% Calculate additional fields from forcing data
if opts.calcFields
    K = length(forcing.t);
    fx = forcing.x;
    fy = forcing.y;
    for k = 1:K
        forcing.div(:,:,k) = divergence(fx,fy,forcing.u(:,:,k),forcing.v(:,:,k));
        [forcing.zeta(:,:,k),~] = curl(fx,fy,forcing.u(:,:,k),forcing.v(:,:,k));
    end
end

end