

forceName = '_input_forcing_anna_long.mat';
forcePath = '../Floe_Cpp/io/inputs/';
fName = strcat(forcePath,forceName);
forcing = load(fName);

%%
K = length(forcing.t);
for k = 1:K
    forcing.div(:,:,k) = divergence( forcing.x, forcing.y, forcing.u(:,:,k), forcing.v(:,:,k) );
    [forcing.zeta(:,:,k),~] = curl( forcing.x, forcing.y, forcing.u(:,:,k), forcing.v(:,:,k) );
    forcing.spd = abs( forcing.u+1i*forcing.v);
end

%%

ftd = forcing.t/86400;
f = 1.4e-4;
ind = 1:10:length(forcing.x);

sqx = 30e3*([0,1]+0.25);
sqy = 30e3*([0,1]+0.5);

figure(1);clf;
colormap(cbrewer2('RdBu',20));

stday = 25;
ks = (24/4)*stday-2;

for k = ks%:K
    pcolor(forcing.x,forcing.y,forcing.spd(:,:,k)); shading flat;
    hold on;
    quiver(forcing.x(ind),forcing.y(ind),forcing.u(ind,ind,k),forcing.v(ind,ind,k),1,'k');
%     plot( sqx([1,1,2,2,1]), sqy([1,2,2,1,1]),'k',LineWidth=1.5)
    hold off;
    caxis(0.25*[-1,1]);
    cb = colorbar;
    daspect([1,1,1]);
    drawnow;
end

%%

x = forcing.x;
y = forcing.y;
t = forcing.t(ks:K); t = t-t(1);
u = forcing.u(:,:,ks:K);
v = forcing.v(:,:,ks:K);

forceName = '_input_forcing_anna_5day.mat';
forcePath = '../Floe_Cpp/io/inputs/';
fName = strcat(forcePath,forceName);
save(fName,'x','y','t','u','v');