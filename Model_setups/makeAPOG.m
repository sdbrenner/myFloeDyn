

% X = 0;
% Y = 0;
% W = 50e3;
% H = 50e3;
% 
% 
% figure(1); clf;
% 
% ax = gca;
% s = 0.35;
% ax.XLim = X+[0+s,1-s]*W;
% ax.YLim = Y+[0+s,1-s]*H;
% 
% 
% text(0.5,0.5,'APOG',...
%      Units='normalized',HorizontalAlignment='center',VerticalAlignment='middle',...
%      FontSize=256, FontWeight='bold', FontName='Century Gothic' );
% daspect([1,1,1]);
% 
% %%
% [Ax,Ay] = ginput;
% [Px,Py] = ginput;
% [Ox,Oy] = ginput;
% [Gx,Gy] = ginput;
% 
% pShape(1) = polyshape( Ax,Ay );
% pShape(2) = polyshape( Px,Py );
% pShape(3) = polyshape( Ox,Oy );
% pShape(4) = polyshape( Gx,Gy );

%%


load('../APOGpoly.mat')
X = 0;
Y = 0;
W = 50e3;
H = 50e3;
window = [X+W*[0,1],Y+H*[0,1]];
% [xc(1),xc(2)] = centroid(pShape(1));
xc = 0.5*[W,H];
for m = 1:length(pShape)
    pShape(m).Vertices = 3*(pShape(m).Vertices-xc ) +xc;
end
pShape(1).Vertices(8:9,1)   = pShape(1).Vertices(8:9,1)-250;
pShape(1).Vertices(end+1,:) = pShape(1).Vertices(1,:)+[250,1];

pShape(2).Vertices(17:19,1) = pShape(2).Vertices(17:19,1)+250;
pShape(2).Vertices(end+1,:) = pShape(2).Vertices(1,:)+[250,1];

pShape(3).Vertices(38:39,1)  = pShape(3).Vertices(38:39,1)+125;
pShape(3).Vertices([1,68],1) = pShape(3).Vertices([1,68],1)-125;
pShape(3).Vertices(end+1,:)  = pShape(3).Vertices(1,:)+[1,500];

pShape(4).Vertices(end+1,:)  = pShape(4).Vertices(1,:)+[250,1];

makeInput(pShape,window,'in_APOG')
% save('../APOGpoly2.mat')


figure(2); clf;

ax = gca;
s = 0.0;
ax.XLim = X+[0+s,1-s]*W;
ax.YLim = Y+[0+s,1-s]*H;

hold on;
for n = 1:4;
%     plot( pShape(n) );
    x = pShape(n).Vertices(:,1);
    y = pShape(n).Vertices(:,2);
%     x = [x;x(end)];
%     y = [y;y(1)+500];
    plot( x, y,'k.-');
    [xc,yc] = centroid(pShape(n));
    plot(xc,yc,'ko');
end
daspect([1,1,1]);