%% ANALYZE STATISTICS OF FLOE SHAPES INCLUDED IN THE FLOEDYN CATALOGUE

%%  Clean workspace and load data
clearvars;
clc
% close all;

load('../../Floe_Cpp/io/inputs/Biblio_Floes.mat');

%%
% G{end+1} = getSquareFloe();
% G{end+1} = getCircularFloe();
N = numel(G);
Atarget = 1;

for n = 1:N
    G{n} = areaNorm( G{n}, Atarget );

    pgon = polyshape( G{n}(:,1), G{n}(:,2) );
%     [xc,yc] = centroid(pgon);
    A(n) = area(pgon);
    P(n) = perimeter(pgon);
    Pch(n) = perimeter( convhull(pgon) );
    [minDiam(n),maxDiam(n),meanDiam(n),cvx(n),thmn(n),thmx(n)] = calliperDiam( pgon );

end
C = 4*pi*A./(P.^2);
[~,ind] = sort(C,'descend');

%% Plot

fH = figure(1); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [8,6]; 

subplot(2,1,1);
plot( 1:N, C(ind),'k');
ylim([0.5,1]);
grid on;
xlabel('Floe #');
ylabel('Circularity');
yyaxis right;
hold on;
plot( 1:N, cvx(ind) );
% plot( 1:N, Pch(ind)./P(ind) );
ylim([0.5,1]);
ylabel('Convexity');


subplot(2,1,2);
hold on;
pH(1) = patch( [1:N,N:-1:1],[minDiam(ind),fliplr(maxDiam(ind))],grey(0.85),EdgeColor='none');
pH(2) = plot( 1:N, meanDiam(ind),'k',LineWidth=2);
pH(3) = plot( 1:N, Pch(ind)./pi,'.',MarkerSize=5 );
% plot( 1:N, sqrt(cvx(ind)).*P(ind)./pi,'b',LineWidth=0.5 )
% plot( 1:N, sqrt(A(ind)) )
% yline( [sqrt(2*Atarget),sqrt(4*Atarget/pi)], '--',Color=grey(0.25),LineWidth=1 )
plot( 1:N, [sqrt(4*A(ind)/pi)], '--',Color=grey(0.25),LineWidth=1 )
% yline( [sqrt(4*Atarget/pi)], '--',Color=grey(0.25),LineWidth=1 )
grid on;
xlabel('Floe #');
ylabel('Floe caliper diameter');

[legH,icons] = legend( pH([1,3,2]), {'Caliper diameter: mean + range','P_{CH}/\pi',''},Location='northwest');
icons(7).YData = ( mean(icons(4).Vertices(1:2,2)) )*[1,1];
icons(6).MarkerSize = 8;
legH.Box = 'off';

%% Plot each floe

F = areaNorm( getCircularFloe, Atarget );
pgonF = polyshape(F);
ln = sqrt(2*Atarget)*[-1,1];

fH = figure(10); clf;
fH.Units = 'inches';
fH.Position([3,4]) = 3;

ax = gca;
hold on;
pgF = plot(pgonF,FaceColor='none',LineStyle='--',LineWidth=1);
pgC = plot(convhull(pgon),FaceColor='none',EdgeColor=grey(0.75),LineWidth=1);
pgG = plot(pgon,LineWidth=1);
lmnH = plot( ln, ln,'-',Color=colourShift( lines(1),+0.35)  );
lmxH = plot( ln, -ln,'-',Color=colourShift( lines(1),-0.35)  );

dataLab = sprintf('min. diam. = %2.2f\nmax. diam. = %2.2f\nC = %2.2f',...
                   minDiam(1),maxDiam(1),C(1) );
tH = text(-0.95*sqrt(Atarget),-0.95*sqrt(Atarget),dataLab);
tH.VerticalAlignment='bottom'; ...
tH.FontName=ax.FontName;
tH.FontSize = 8;
xlim(sqrt(Atarget)*[-1,1]);
ylim(sqrt(Atarget)*[-1,1]);
daspect([1,1,1]);
imFrame = cell(1,N);

for n = 1:N
    m = ind(n);
    pgon = polyshape(G{m});
    pgG.Shape = pgon;
    pgC.Shape = convhull(pgon);
    lmnH.XData  = ln*cos(thmn(m));
    lmxH.XData  = ln*cos(thmx(m));
    lmnH.YData  = ln*sin(thmn(m));
    lmxH.YData  = ln*sin(thmx(m));
    dataLab = sprintf('catalog num. = %g\nmin. diam. = %2.2f,  max. diam. = %2.2f\ncircularity = %2.2f,  convexity = %2.2f',...
                       m, minDiam(m),maxDiam(m),C(m),cvx(m) );
    tH.String = dataLab;
    drawnow;
    imFrame{n} = getframe;
    pause(0.2);
%     waitforbuttonpress;
end

%% Test
% 
% n = 36;
% pgon = polyshape( G{n} );
% 
% N = 1e3;
% th = linspace(0,2*pi,N);
% [xc,yc] = centroid(pgon);
% V = pgon.Vertices-[xc,yc];
% D = zeros(1,N);
% for n = 1:N
%     R = [ cos(th(n)), -sin(th(n)) ;
%           sin(th(n)),  cos(th(n))  ];
%     VV = [V(:,1).';V(:,2).']';
%     VR = VV*R;
%     D(n) = range(VR(:,1));
%     Pn(n) = convexity( VR );
% end
% [minD,ind] = min(D);   thetaMinD = th(ind);
% [maxD,ind] = max(D); thetaMaxD = th(ind);
% 
% figure;
% plot( th, D );
% xline(thetaMinD);
% xline(thetaMaxD);

%%


function [minD,maxD,meanD,P,thetaMinD,thetaMaxD] = calliperDiam( pgon )
    N = 1e3;
    th = linspace(0,pi,N);
    [xc,yc] = centroid(pgon);
    V = pgon.Vertices-[xc,yc];
    D = zeros(1,N);
    for n = 1:N
        R = [ cos(th(n)), -sin(th(n)) ;
              sin(th(n)),  cos(th(n))  ];
        VV = [V(:,1).';V(:,2).']';
        VR = VV*R;
        D(n) = range(VR(:,1));
        Pn(n) = convexity( VR );
    end
    [minD,ind] = min(D);   thetaMinD = th(ind);
    [maxD,ind] = max(D); thetaMaxD = th(ind);
    meanD = mean(D);
    P = min(Pn);
end

function  floeXY = getCircularFloe()
    theta = linspace(0,2*pi,360);
    [x,y] = pol2cart(theta,100);
    floeXY = [x(:),y(:)];
end

function  floeXY = getSquareFloe()
    D = 100;
    x = [1,1,2,2];
    y = [1,2,2,1];
    floeXY = D*[x(:),y(:)];
end

function F = areaNorm( G, Atarget )

    pgon = polyshape( G(:,1), G(:,2) );
    A = area(pgon);
    [xc,yc] = centroid(pgon);
    F = (G-[xc,yc]) * sqrt(Atarget/A);
%     F = (G-[xc,yc]) * (pi*Atarget/perimeter(convhull(pgon))); % MCD target
    
end



function P = convexity( VR )
    % Definition per Zunic & Rosin, 2002
    VR(end+1,:) = VR(1,:);
    P2 = 2*sum(range(VR));
    P1 = sum(sum(abs(diff(VR))));
    P = P2/P1;
end
