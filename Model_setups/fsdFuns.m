function [fx,Fx] = fsdFuns( minSize,maxSize,alpha ) 
% FSDFUNS generates anonymous floe-size-distribution functions for bounded
% power law FSDs (following Stern et. al., 2018)
%
%   [fx,Fx] = fsdFuns( minSize,maxSize,alpha ) 
%   gives:
%       fx the non-cumulative distrubution function (probability density function); and,
%       Fx the complementary cumulative distribution function
%
%   For more details, see:
%       On reconciling disparate studies of the sea-ice floe size distribution
%       Stern, Schweiger, Zhang, and Steele, 2018
%       https://doi.org/10.1525/elementa.304
%
%   S.D.Brenner, 2022

    % parameters associated with distribution 
    B = (maxSize/minSize).^(1-alpha);
    c = (alpha-1)*minSize^(alpha-1)/(1-B);
    C = c/(alpha-1);
    R = C*maxSize^(1-alpha);

    % Theoretical FSD functions (see Stern2018)
    fx = @(x) c*x.^(-alpha);
    Fx = @(x) C*x.^(1-alpha)-R;

end