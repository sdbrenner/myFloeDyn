function [floeOut,forcing] = calcIOBLDiss( floeOut, forcing )

%   [floeOut,forcing] = calcIOBLDiss( floeOut, forcing )
%
%   Note: code vectorization is possible instead of the nested for-loops,
%   but some experimentation has shown that the large matrices needed for
%   thatm approach lead to even slower code execution (at least for me).
%   The current approach is still not ideal and I feel there must be some
%   more clever steps I can take to generate the full ice velocity field.

    %% Create floe masks if necessary

    if ~isfield( floeOut.floes,'floeMasks')
        floeOut = makeFloeMask( floeOut, forcing.x,forcing.y );
    end

    %% Generate full ice velocity field
    
    [M,L] = size(floeOut.floes.floeMasks);
    [X,Y] = meshgrid( forcing.x, forcing.y );
    W = range( floeOut.chunk.win(1:2));
    H = range( floeOut.chunk.win(3:4));
    pbcXY = @(X,xc,W) (X-xc) + W*((xc-X)>0.5*W) -W*((X-xc)>0.5*W);
        
    nx = length(forcing.x); ny = length(forcing.y);
    Ui = NaN( [ny,nx,L] );
    disp('Calculating full ice velocity field...')
    % Loop through time
    for l = 1:L
        fprintf('timestep %2g of %2g...\n',l,L);
        UiL = NaN( [ny,nx] );
        % Loop through floes
        for m = 1:M
            % translational ice velocity
            UiT = floeOut.state.ui(m,l) + 1i*floeOut.state.vi(m,l);
            iceMaskML = logical(floeOut.floes.floeMasks{m,l});
            % rotational ice velocity
            xc = floeOut.state.xp(m,l);
            yc = floeOut.state.yp(m,l);
            omega = floeOut.state.om(m,l);
%             xr = (forcing.x-xc) + W*((xc-forcing.x)>0.5*W) -W*((forcing.x-xc)>0.5*W);
%             yr = (forcing.y-yc) + H*((yc-forcing.y)>0.5*H) -H*((forcing.y-yc)>0.5*H);            
%             UiR = omega*( -yr.' + 1i*xr );
            UiR = (1i)*iceMaskML; % allocate
            UiR(iceMaskML) = omega*( -pbcXY(Y(iceMaskML),yc,H) +1i*pbcXY(X(iceMaskML),xc,W) );
            % Add total ice velocity for floe to full matrix
            UiL(iceMaskML) = UiT + UiR(iceMaskML);       
        end
        Ui(:,:,l) = UiL;
    end
    
    
    %% Calculate dissipation
    
    % Extract forcing data
    Uo = forcing.u(:,:,1:L) + 1i*forcing.v(:,:,1:L);
    if isfield( forcing, 'ug')
        Ug = forcing.ug(:,:,1:L) + 1i*forcing.vg(:,:,1:L);
    else
        Ug = Uo;
    end

    % Define constant parameters
    rho = 1025;
    Cio = 5e-3;
    f = 1.4e-4;
    d = 1;
    
    % Calculate stress/dissipation
    Urel = Uo-Ui;
    Tio = rho*Cio*Urel.*abs(Urel);
    diss = real(Tio.*conj(Urel)); 
    
    % Calculate equivalent steady-state dissipation
    UiSS = seaIceSteadyState( Cio/(f*d), Uo, Ug );
    steadyDiss = rho*Cio*abs(UiSS-Uo).^3;
    
    %% Integrate and Average over floes
    
    dx = mean(diff(forcing.x));
    dy = mean(diff(forcing.y));
    % Pre-allocate memory
    floeOut.chunk.ioblDiss = NaN(1,L);
    floeOut.ocean.avg.diss = NaN(M,L);
    floeOut.ocean.avg.excessDiss = NaN(M,L);

    for l  = 1:L
        dissL = diss(:,:,l);
        steadyDissL = steadyDiss(:,:,l);

        % total integrated IOBL dissipation
        floeOut.chunk.ioblDiss(l) = dx*dy*sum(dissL,'all','omitnan');
        floeOut.chunk.ioblSteadyDiss(l) = dx*dy*full(sum(floeOut.chunk.iceMask{l}.*steadyDissL,'all','omitnan'));
        % Loop through floes
        for m = 1:M
            iceMaskML = floeOut.floes.floeMasks{m,l};
            
            floeOut.ocean.avg.diss(m,l) = mean( dissL(iceMaskML) );
            floeOut.ocean.avg.steadyDiss(m,l) = mean( steadyDissL(iceMaskML) );
            floeOut.ocean.avg.excessDiss(m,l) = mean( dissL(iceMaskML)-steadyDissL(iceMaskML) );
            floeOut.ocean.avg.excessDissRat(m,l) = mean( dissL(iceMaskML)./steadyDissL(iceMaskML) );
        end
        
    end

    %% Add fields to forcing file

    forcing.ui = real(Ui);
    forcing.vi = imag(Ui);
    forcing.diss = diss;
    forcing.steadyDiss = steadyDiss;


end
