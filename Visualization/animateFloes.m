function imFrame = animateFloes( floeOut,forcing,figOpts )



%% Parse input arguments
% check for all required fields in floeOut
% check if 'forcing' exists or not
% check for figOpts (and assign defaults)


% Default figure options
arguments
    floeOut             {mustBeA(floeOut,"struct")}
    forcing             %{mustBeA(forcing,"struct")}
    figOpts.figNum      = 2;
    figOpts.fontName    = 'Open Sans'; %(SASIP font)
    figOpts.fontSize    = 18;
    figOpts.resolution  = 720;
    figOpts.aspectRatio = 1;
    figOpts.iceColour   = [1,1,0.95]; %0.95*[1,1,1]; %[0.95,0.95,0.9];
    figOpts.iceAlpha    = 1; 
    figOpts.iceBorder   = 'none'; 
    figOpts.cmap        = flipud(colorcet('L17','N',20));
    figOpts.pFun        = @(X,Y,U,V) abs(U+1i*V);  % speed
    figOpts.CLim        = 0.2*[0,1];               % speed
    figOpts.cbLab       = 'ocean surface current speed [m/s]';
    figOpts.axColour    = [0,0.16539,0.65755];
end

% % overwrite figure options from input
% if nargin>2 && isstruct(figOptsIn)
%     flds = fields(figOptsIn);
%     for n = 1:length(flds)
%         figOpts.(flds{n}) = figOptsIn.(flds{n});
%     end
% end



%% Initialize figure

    
    fH = figure(figOpts.figNum); clf;
    fH.Color = 'w';
    fH.Position([3,4]) = figOpts.resolution*[ figOpts.aspectRatio,1 ];

    ax = gca;
    hold on;
    k = 1;
    
    if ~isempty(forcing)
        forceHandles = plotForcing(forcing,figOpts);
    end
    pgH = plot( floeOut.chunk.icePoly{k} );
    pgH.FaceColor = figOpts.iceColour;
    pgH.FaceAlpha = figOpts.iceAlpha;
    pgH.EdgeColor = figOpts.iceBorder;
    pgH.LineWidth = 1;
    
    
    % Set up axis
    ax.Box = 'on';
    ax.Layer = 'top';
    ax.LineWidth = 1;
    ax.XLim = floeOut.chunk.win(1:2);
    ax.YLim = floeOut.chunk.win(3:4);
    ax.CLim = figOpts.CLim;

    daspect([1,1,1]); 
    
    ax.XAxis.Exponent = 0;
    ax.XAxis.TickLabelFormat = '%.0f';
    ax.YAxis.Exponent = 0;
    ax.YAxis.TickLabelFormat = '%.0f';
    
    ax.FontSize = figOpts.fontSize;
    ax.FontName = figOpts.fontName;
    ax.LineWidth = 1.5;
    XT = linspace( ax.XLim(1),ax.XLim(2),6 );
    XTL = string(XT/1000);
    ax.XTick = XT; 
    ax.XTickLabel = XTL;
    % XT = ax.XTick;
    % XTL = ax.XTickLabel;
    ax.YTick = XT;
    ax.YTickLabel = XTL;
    xlabel('{\itx}-distance [km]'); ylabel('{\ity}-distance [km]')
    
    
    ax.Color = figOpts.axColour;
    colormap( figOpts.cmap );
    
    
    % blurb = {'FloeDyn simulation of ice floe response to submesoscale eddies';...
    %          '(one-way coupled)';...
    blurb = ['Samuel Brenner, ',datestr(now,'yyyy')];
    tH = text(1/50,1/50,blurb,...
              Units="normalized",...
              HorizontalAlignment='left',...
              VerticalAlignment='bottom');
    tH.Color = colourShift( figOpts.cmap(1,:), 0.35 );
    tH.FontAngle='italic';
    tH.FontName = figOpts.fontName;
    tH.FontSize = (2/3)*figOpts.fontSize;
    % try to move text below floes?
    ax.Children([1,2]) = ax.Children([2,1]);
    
    %% Animate
    
    K = length(floeOut.chunk.icePoly);
    imFrame = cell(1,K);
    
    for k = 1:K
        pgH.Shape = floeOut.chunk.icePoly{k};
        if ~isempty(forcing)
            updateForcing(k,forceHandles,forcing,figOpts)
        end        
        tLab = sprintf('time = %sd-%s',...
                       datestr(seconds(floeOut.state.time(k)),"dd" ),...
                       datestr(seconds(floeOut.state.time(k)),"HH:MM" ));
        title( tLab,FontWeight="normal" );
        drawnow;
        imFrame{k} = getframe(fH); 
    end 



end

%%


function forceHandles = plotForcing(forcing,figOpts)
    % Extract data
    k = 1;
    uk = forcing.u(:,:,k);
    vk = forcing.v(:,:,k);
    cdata = figOpts.pFun( forcing.x, forcing.y, uk, vk );
    % Background pcolor
    forceHandles.pcH = pcolor(forcing.x,forcing.y,cdata); shading interp;
    % Quiver
    qInd = 1:15:length(forcing.x);
    arSc = 5000;
    forceHandles.qH = quiver( forcing.x(qInd),forcing.y(qInd),...
                              arSc*uk(qInd,qInd), arSc*vk(qInd,qInd), 0 );
    forceHandles.qH.Color = colourShift( figOpts.cmap(6,:), -0.35 );
    forceHandles.qH.MaxHeadSize = 1;
    % Make colourbar
    cb = colorbar;
    cb.LineWidth = 1;
    cb.FontName = figOpts.fontName;
    ylabel(cb,figOpts.cbLab);

end

function updateForcing(k,forceHandles,forcing,figOpts)
    % Extract data
    uk = forcing.u(:,:,k);
    vk = forcing.v(:,:,k);
    cdata = figOpts.pFun( forcing.x, forcing.y, uk, vk );
    
    % Background pcolor
    forceHandles.pcH.CData = cdata;
%     pcH.CData = zetak/f; caxis(3.1*[-1,1]); ylabel(cb,'ocean surface vorticity (\zeta/f)');
%     pcH.CData = divk; caxis(0.5*3.1e-4*[-1,1]); ylabel(cb,'ocean surface divergence [1/s]');
    % Quiver
    qInd = 1:15:length(forcing.x);
    arSc = 5000;
    forceHandles.qH.UData = arSc*uk(qInd,qInd); 
    forceHandles.qH.VData = arSc*vk(qInd,qInd);
end






