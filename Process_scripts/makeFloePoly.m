function floeOut = makeFloePoly( floeOut )
% MAKEFLOEPOLY returns the coordinates MATLAB polyshapes corresponding to
% each of the floes
%
%   floeOut = makeFloePoly( floeOut ), where 'floeOut' was produced by the
%   'readFloeOut' function. The 'floeOut' structure is modified to include
%   a new field:
%       floeOut.floes.floePoly
%
%   S.D.Brenner, 2022


% Check if the floeXYT field already exists; if not, run calcFloeXYT
if ~ isfield( floeOut.floes, 'floeXYT' )
    floeOut = calcFloeXYT( floeOut );
end