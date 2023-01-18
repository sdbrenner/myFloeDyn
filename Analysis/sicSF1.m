
clear; 
clc;
% close all;

batchCase = 'FSD_300_2500_2_geo';
fName = strcat('floeOut_',batchCase );
load( strcat('../MAT_Files/',fName) );
n = 2;
%%
% Generate polyshape of intial and final sea ice configuration
L = length( floeOut(n).state.time );
fol = subsetFloeOut( floeOut(n), [1,L] );
fol = makeDomainPoly( fol );


%% Randomly generate points
X = fol.chunk.win(1);
W = range(fol.chunk.win(1:2));
Y = fol.chunk.win(3);
H = range(fol.chunk.win(3:4));

rng(42);
numPts = 5000;
ptx = X + rand(1,numPts)*W;
pty = Y + rand(1,numPts)*H;
meanPtSpacing = sqrt( W*H/numPts );
dists = ToroidalDistance(ptx,pty,fol.chunk.win);

% Assign status
for l = 1:2
    % sic = isinterior( fol.chunk.icePoly{l},ptx,pty );
    sic = isIn( fol, ptx,pty, l );
    
    % Calculate distances and SF1
    % dx = (ptx-ptx');
    % dy = (ptx-ptx');
    % dists = sqrt( dx.^2 + dy.^2 );
%     dists = ToroidalDistance(ptx,pty,fol.chunk.win);
    SF1 = abs(sic-sic.');
    SF1i = SF1 + 1i*(sic & sic');
%     SF1ii = SF1 -(sic & sic');
%     ic = double(sic); ic(~ic)=1i;
%     SF1i = 1+( real(ic)-imag(ic.') );
%     SF1i(SF1i==2) = 1i;
    
    % Bin-average
    binSize = 200;
    bins = 0:binSize:W/2; %( sqrt(W^2+H^2)/2 );
%     bins = logspace( log10( 100 ), log10( mean([W,H]) ), 20 );
    [~,bSF1] = binFun( dists(:), SF1(:) ,bins,'returnX','binCenters' );
    [bDist,bSF1i] = binFun( dists(:), SF1i(:) ,bins,'returnX','binCenters' );
%     [bDist,bSF1ii] = binFun( dists(:), SF1ii(:) ,bins,'returnX','binCenters' );
    
    %% Plot
    fH = figure(1); if l==1, clf; end
    fH.Units = 'inches';
    fH.Position([3,4])=6;
%     colormap( colorcet('CBD2') )
    cols = lines(2);

    ax = subplot(2,2,l); cla;
    scatter( ptx, pty, 2, sic,'filled',MarkerFaceAlpha=0.35);
    hold on;
    pgH = plot(fol.chunk.icePoly{l},FaceColor='none',EdgeColor=grey(0.35),LineWidth=0.5 );
    hold off;
    xlim(X+[0,W])
    ylim(Y+[0,H]);
    daspect([1,1,1]);
    colormap( ax, [cols(l,:);grey(0.65)] )
    tLab = ["t = 0";"t = 72 h"];
    title(tLab(l));

    ax = subplot(2,2,3:4); 
    hold on;
    SIC = mean(sic(:));
%     plot( bDist, 2*(SIC-imag(bSF1i)), Color=grey(0.85),LineWidth=2.5)
    plot( bDist, real(bSF1i),'-',LineWidth=1.5,Color=cols(l,:) );
%     plot( bDist, imag(bSF1i),'--',LineWidth=1.5,Color=cols(l,:) );
%     plot( bDist, real(bSF1ii),'-',LineWidth=1.5,Color=cols(l,:) );
%     plot( bDist, 2*SIC-3*imag(bSF1i),'--',Color=cols(l,:) );

%     plot( bDist, 0.5*real(bSF1i)+imag(bSF1i),'-.',LineWidth=1.5,Color=cols(l,:) );
    ax.XScale = 'log'; ax.YScale = 'log';
    ax.XLim = [100,W/2];%sqrt(W^2+H^2)/2];
    ax.YLim = [3e-2,1e0];
%     grid on;
    ylabel('first-order structure function');
    xlabel('distances [m]');
    if l==2  
%         yline( [ SIC, SIC^2, 2*SIC*(1-SIC) ],'-',...
%                Color=grey(0.65),LineWidth=1); 
        yline( [ -SIC, 2*SIC-3*SIC.^2 ],'-',...
               Color=grey(0.65),LineWidth=1); 
        xline( max(fol.floes.floeSize),'k--',Color=grey(0.35),LineWidth=1); 
        xline( meanPtSpacing,'k:' );
    end
%     ax.XScale = 'linear'; 
    ax.YScale = 'linear'; ax.YLim = [0,1]; ax.YLim = [-1,1];

end

%% hit test
function in = isIn(fol,ptx,pty,l)
    
    % [ BASIC ]
    M = fol.floes.numFloes;
    pbcheck = zeros(1,M);
    inm = zeros( length(ptx), M );
    for m = 1:M
        FX = fol.floes.floeXYT{m,l}(1,:);
        FY = fol.floes.floeXYT{m,l}(2,:);
        inm(:,m) = inpolygon( ptx, pty, FX,FY);
        % check if polygon overlaps boundary
        pbcheck(m) = any( ~isbetween(FX,fol.chunk.win(1:2)) |...
                          ~isbetween(FY,fol.chunk.win(3:4)) );
    end
    % [ ACCOUNTING FOR PBC ]
    J = sum(pbcheck); % number of floes that overlap boundaries
    inj = zeros( length(ptx), 9*J );
    pbind = find(pbcheck);
    W = range(fol.chunk.win(1:2)); H = range(fol.chunk.win(3:4));
    xOff = W*reshape(repmat(-1:1,3,1) ,1,[]);
    yOff = H*reshape(repmat(-1:1,3,1)',1,[]);
    k = 1;
    for j = 1:J
        m = pbind(j);
        for r = 1:9
            FX = fol.floes.floeXYT{m,l}(1,:) + xOff(r);
            FY = fol.floes.floeXYT{m,l}(2,:) + yOff(r);
            inj(:,k) = inpolygon( ptx, pty, FX,FY);
            k = k+1;    
        end
    end


    in = any([inm,inj],2);
end


function dists = ToroidalDistance(ptx,pty,win)
    
    X = win(1); Y = win(3);
    W = range(win(1:2)); H = range(win(3:4));

    % transform from X+[0,W]->[0,W]
    ptx = ptx-X;
    pty = pty-Y;
    
    % Calculate regular dx dy for each point set
    dx = abs(ptx-ptx');
    dy = abs(pty-pty');
    
    % Modify dx dy
    dx(dx>0.5*W) = W-dx(dx>0.5*W);
    dy(dy>0.5*H) = H-dy(dy>0.5*H);
    
    % now calculate distance
    dists = sqrt( dx.^2 + dy.^2 );

end






% float ToroidalDistance (float x1, float y1, float x2, float y2)
% {
%     float dx = std::abs(x2 - x1);
%     float dy = std::abs(y2 - y1);
%  
%     if (dx > 0.5f)
%         dx = 1.0f - dx;
%  
%     if (dy > 0.5f)
%         dy = 1.0f - dy;
%  
%     return std::sqrt(dx*dx + dy*dy);
% }

