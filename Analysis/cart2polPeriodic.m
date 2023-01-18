function [TH,R] = cart2polPeriodic(X,Y,win)
% CAR2POLPERIODIC Transform Cartesian coordinates in a periodic domain to
% polar coordinates.
%
%   [TH,R] = car2polPeriodic(X,Y,win)  transforms corresponding elements of
%   data stored in Cartesian coordinates X,Y in the periodic domain, 'win'
%   (x1,x2,y1,y2), to polar coordinates (angle TH and radius R).
%
%   S.D.Brenner, 2022

    arguments 
        X   {mustBeNumeric}
        Y   {mustBeNumeric}
        win (1,4) {mustBeNumeric} = [min(X(:)),range(X(:)),min(Y(:)),range(Y(:))]
    end

    W = range( win(1:2) );
    H = range( win(3:4) );
    Xr = X + W*((-X)>0.5*W) -W*((X)>=0.5*W);
    Yr = Y + H*((-Y)>0.5*H) -H*((Y)>=0.5*H);
    R = hypot( Xr,Yr );
    TH = atan2( Yr,Xr ); 

end


%% DEPRICATED VERSION (slower and more complicated)
%{
    % Get domain information
%     X0 = win(1);
    W = range( win(1:2) );
%     Y0 = win(3);
    H = range( win(3:4) );
    
    % Coordinates rescale functions (convert to angles on circular domains)
%     th = @(x) (x-X0)*(2*pi)/W -pi;
%     ph = @(y) (y-Y0)*(2*pi)/H -pi;
    
    % Get angle differences
%     dth = atan2( sin(th(X)-th(0)), cos(th(X)-th(0)) );
%     dph = atan2( sin(ph(Y)-ph(0)), cos(ph(Y)-ph(0)) );
    dth = atan2( sin(2*pi*X/W), cos(2*pi*X/W) );
    dph = atan2( sin(2*pi*Y/H), cos(2*pi*Y/H) );


    % Convert to distances (i.e., arclengths on circular domains)
    dx = (dth)*W/(2*pi);
    dy = (dph)*H/(2*pi);
    
    % Calculate global radial distance and angle
    % [NOTE: I don't know if this is the right way to define theta -- it leads
    % to discontinuities at the wraparound points. I need to investigate this
    % further ]
    R = hypot( dx,dy );
    TH = atan2( dy,dx ); 

end
%}
