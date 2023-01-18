function [ui] = seaIceSteadyState(Cstar,uo,ug)
% SEAICESTEADYSTATE calculates the steady-state sea ice velocity in
% response to ocean currents (no wind)
%
%   ui = seaIceSteadyState(Cstar,uo)
%   ui = seaIceSteadyState(Cstar,uo,ug) where: Cstar = Cio/(di*f) [scalar],
%   uo is the ocean velocity, and ug is the (optional) geostrophic
%   velocity, solves the steady-state sea ice momentum equation subject to
%   Coriolis accleration, ice-ocean stress, and sea surface tilt (expressed
%   as a geostrophic velocity):
%
%       rhoo*di*i*f*(ui-ug) = rhoo*Cio*(ui-uo)*||ui-uo||
%
%   S.D.Brenner, 2022


    if nargin<3 || isempty(ug), ug=0; end
    
    % Ageostrophic velocity:
    uag = uo-ug; 
    
    % Solve for turning angle
    Y = Cstar*abs(uag);
    beta = 2*atan(sqrt(sqrt(4*Y.^2+1)-2*Y));

    % Calculate velocity
    ui = ug + uag.*cos(beta).*exp(-1i*beta);

    %{ 
    % [ Depricated: I figured out an analytical solution based on geometic
    arguments ]
    %
    % if the inputs are "large", don't compute the ice velocities directly;
    % instead generate an interpolation vector
    N = 250;
    if numel(uag)>N 
        % interpolator:
%         uos = linspace(0,0.3,N);
        uos = [0,logspace(-4,-0.3,N-1)];
        uis = ssIceFun( Cstar,uos );  
        uig = interp1( uos, uis, abs(uag),'pchip','extrap' ); % uig in uag frame-of-reference
        % rotatation back to regular frame of reference
        uig = uig.*exp(1i*angle(uag));
    else
        uig = ssIceFun( Cstar,uag );
    end

    ui = uig+ug;
    err = abs( 1i*(ui-ug)-(Cstar*(uo-ui).*abs(uo-ui)) );
    %}


end


%     if nargin<3 || isempty(ug), ug = 0*uo; end
%     if length(ug)~=1, ug = repmat(ug,size(uo)); end
%     
%     ui = NaN(size(uo));
%     err = ui;
%     syms uin;
%     
%     for n = 1:numel(uo)
%         eqn = 1i*(uin-ug(n))==Cstar*(uin-uo(n)).*abs(uin-uo(n));
%         try
%             ui(n) = vpasolve(eqn,uin,uo(n));
%             err(n) = subs(eqn,ui(n));
%         catch
%         end
%     end
% 



function uis = ssIceFun( Cstar,uos )
    % normalized units (t-star ~ t*f)
    % NB: there's almost certainly a smarter/faster way to do this (e.g.,
    % Newton–Raphson), but this seems to work okay
    
    dt = 1e-3;    
    M = numel(uos);
    uis = zeros( size(uos) );
    for m = 1:M
        uo = uos(m);
        T = min([20,2/abs(uo)]);
        t = 0:dt:T;
        N = length(t);
        ui = zeros(1,N);
        for n = 1:N-1
            Tio = Cstar*(ui(n)-uo)*abs(ui(n)-uo);
            ui(n+1) = ui(n) -dt*( Tio + 1i*ui(n)  );
        end
        uis(m) = ui(end);
    end
end







