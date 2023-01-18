%% P-S Separation
% Following Smith, 2007

%% clean workspace

clearvars -except forcing;
clc;
% close all;

%% Load forcing

forceName = '_input_forcing_SM_5day';
forcePath = '../Floe_Cpp/io/inputs/forcings/';
fName = strcat(forcePath,forceName);
if ~exist('forcing','var')
    forcing = load(fName);
end


%% Loop through time and perform decomposition

[N,M,L] = size(forcing.u);
x = forcing.x;
y = forcing.y;

forcingP = forcing; 
forcingS = forcing;

for l = 1:L
    % Extract data from structure
    u = forcing.u(:,:,l);
    v = forcing.v(:,:,l);
    [zeta,~] = curl(x,y,u,v);
    [div] = divergence(x,y,u,v);
    % decompose
    [Pu,Pv,Su,Sv] = psdecomp(u,v);

    % add to new forcing files
    forcingP.u(:,:,l) = Pu;
    forcingP.v(:,:,l) = Pv;
    forcingS.u(:,:,l) = Su;
    forcingS.v(:,:,l) = Sv;
end

%% Plot animations


% loop through time
for l = L
    % Extract data from structure
    u = forcing.u(:,:,l);
    v = forcing.v(:,:,l);
    [zeta,~] = curl(x,y,u,v);
    [div] = divergence(x,y,u,v);
    % Potential
    Pu = forcingP.u(:,:,l);
    Pv = forcingP.v(:,:,l);
    [zetaP,~] = curl(x,y,Pu,Pv);
    [divP] = divergence(x,y,Pu,Pv);    
    % Solenoidal
    Su = forcingS.u(:,:,l);
    Sv = forcingS.v(:,:,l);
    [zetaS,~] = curl(x,y,Su,Sv);
    [divS] = divergence(x,y,Su,Sv);
    % Reconstitute (to check);
    Ru = Pu+Su;
    Rv = Pv+Sv;
    [zetaR,~] = curl(x,y,Ru,Rv);
    [divR] = divergence(x,y,Ru,Rv);


    % plot parameters
    CL = 1.4e-4*[-1,1];
    XL = [min(x),max(x)];
    YL = [min(y),max(y)];
    ssInd = 1:20:length(x);

    % plot vorticity
    fH = figure(1); clf;
    fH.Units = 'inches'; fH.Position([3,4]) = [10,3];
    subplot(1,3,1);
    hold on;
    pcolor(x,y,zeta); shading interp;
    quiver(x(ssInd),y(ssInd),u(ssInd,ssInd),v(ssInd,ssInd),1,'k');
    caxis(CL);
    xlim(XL);
    ylim(YL);
    colormap(cbrewer2('RdBu',21));
    title('Full');

    subplot(1,3,2);
    hold on;
    pcolor(x,y,zetaP); shading interp;
    quiver(x(ssInd),y(ssInd),Pu(ssInd,ssInd),Pv(ssInd,ssInd),1,'k');
    caxis(CL);
    xlim(XL);
    ylim(YL);    
    colormap(cbrewer2('RdBu',21));
    title('Potential')

    subplot(1,3,3);
    hold on;
    pcolor(x,y,zetaS); shading interp;
    quiver(x(ssInd),y(ssInd),Su(ssInd,ssInd),Sv(ssInd,ssInd),1,'k');
    caxis(CL);
    xlim(XL);
    ylim(YL);
    colormap(cbrewer2('RdBu',21));
    title('Solenoidal')

%     subplot(2,2,4);
%     hold on;
%     pcolor(x,y,zetaR); shading interp;
%     quiver(x(ssInd),y(ssInd),Ru(ssInd,ssInd),Rv(ssInd,ssInd),1,'k');
%     caxis(CL);
%     xlim(XL);
%     ylim(YL);
%     colormap(cbrewer2('RdBu'));
%     title('Reconstruction');
    sgtitle('Vorticity');
%     drawnow;

    % plot divergence
    fH = figure(2); clf;
    fH.Units = 'inches'; fH.Position([3,4]) = [10,3];
    subplot(1,3,1);
    hold on;
    pcolor(x,y,div); shading interp;
    quiver(x(ssInd),y(ssInd),u(ssInd,ssInd),v(ssInd,ssInd),1,'k');
    caxis(CL);
    xlim(XL);
    ylim(YL);
    colormap(cbrewer2('RdBu',21));
    title('Full');

    subplot(1,3,2);
    hold on;
    pcolor(x,y,divP); shading interp;
    quiver(x(ssInd),y(ssInd),Pu(ssInd,ssInd),Pv(ssInd,ssInd),1,'k');
    caxis(CL);
    xlim(XL);
    ylim(YL);    
    colormap(cbrewer2('RdBu',21));
    title('Potential')

    subplot(1,3,3);
    hold on;
    pcolor(x,y,divS); shading interp;
    quiver(x(ssInd),y(ssInd),Su(ssInd,ssInd),Sv(ssInd,ssInd),1,'k');
    caxis(CL);
    xlim(XL);
    ylim(YL);
    colormap(cbrewer2('RdBu',21));
    title('Solenoidal')

%     subplot(2,2,4);
%     hold on;
%     pcolor(x,y,divR); shading interp;
%     quiver(x(ssInd),y(ssInd),Ru(ssInd,ssInd),Rv(ssInd,ssInd),1,'k');
%     caxis(CL);
%     xlim(XL);
%     ylim(YL);
%     colormap(cbrewer2('RdBu'));
%     title('Reconstruction');
    sgtitle('Divergence');
%     drawnow;
end



%% Plot velocities

% loop through time
for l = 20%:L
    % Extract data from structure
    u = forcing.u(:,:,l);
    v = forcing.v(:,:,l);
    [zeta,~] = curl(x,y,u,v);
    [div] = divergence(x,y,u,v);
    % Potential
    Pu = forcingP.u(:,:,l);
    Pv = forcingP.v(:,:,l);
    [zetaP,~] = curl(x,y,Pu,Pv);
    [divP] = divergence(x,y,Pu,Pv);    
    % Solenoidal
    Su = forcingS.u(:,:,l);
    Sv = forcingS.v(:,:,l);
    [zetaS,~] = curl(x,y,Su,Sv);
    [divS] = divergence(x,y,Su,Sv);
    % Reconstitute (to check);
    Ru = Pu+Su;
    Rv = Pv+Sv;
    [zetaR,~] = curl(x,y,Ru,Rv);
    [divR] = divergence(x,y,Ru,Rv);


    % plot parameters
    XL = [min(x),max(x)]/1000;
    YL = [min(y),max(y)]/1000;
    ssInd = 1:30:length(x);

    fH = figure(3); clf;
    fH.Units = 'inches';
    fH.Position([3,4]) = [10,3];
    subplot(1,3,1);
    hold on;
    quiver(x(ssInd)/1000,y(ssInd)/1000,u(ssInd,ssInd),v(ssInd,ssInd),1,'k');
    caxis(CL);
    xlim(XL);
    ylim(YL);
    colormap(cbrewer2('RdBu'));
    title('Full');
    daspect([1,1,1]);

    subplot(1,3,2);
    hold on;
    quiver(x(ssInd)/1000,y(ssInd)/1000,Pu(ssInd,ssInd),Pv(ssInd,ssInd),1,'k');
    caxis(CL);
    xlim(XL);
    ylim(YL);    
    colormap(cbrewer2('RdBu'));
    title('Potential (irrotational)')
    daspect([1,1,1]);

    subplot(1,3,3);
    hold on;
    quiver(x(ssInd)/1000,y(ssInd)/1000,Su(ssInd,ssInd),Sv(ssInd,ssInd),1,'k');
    caxis(CL);
    xlim(XL);
    ylim(YL);
    colormap(cbrewer2('RdBu'));
    title('Solenoidal (nondivergent)')
    daspect([1,1,1]);

end

%% Save


fName = strcat(forcePath,forceName);
% save(strcat(fName,'_irrot'),'-struct','forcingP')
% save(strcat(fName,'_nondiv'),'-struct','forcingS')

%% Function

function [Pu,Pv,Su,Sv] = psdecomp(uu,vv)
    % Separate (uu,vv) into irrotational (potential flow) and 
    % nondivergent (Solenoidal flow) parts. This implementation 
    % assumes dx = dy, though the arrays can be rectangular. 
    % original implementation by J. A. Smith, May 23, 2005; 
    % function updated by J. A. Smith, June 28, 2007. 
    
    [mn, nn] = size(uu); % vv is the same size. 
%     nft = 2^nextpow2(max(mn,nn)*1.2);% min 20% zero pad. 
    nft = max(mn,nn); % don't zero pad
    nf2 = nft/2; 
    kv = (-nf2:(nf2-1))/nft; 
    thetamn= fftshift(atan2(kv(ones(nft,1),:)',kv(ones(nft,1),:))); 
    vk1 = fft2(uu,nft,nft); 
    vk2 = fft2(vv,nft,nft); 
    vkdiv = vk1.*cos(thetamn)+vk2.*sin(thetamn); 
    vav = real(ifft2(vkdiv.*cos(thetamn))); 
    Pu = vav(1:mn,1:nn); 
    vav = real(ifft2(vkdiv.*sin(thetamn))); 
    Pv = vav(1:mn,1:nn); 
    vk1 = fft2(uu,nft,nft); 
    vk2 = fft2(vv,nft,nft); 
    vkdiv = vk2.*cos(thetamn)-vk1.*sin(thetamn); 
    vav = real(ifft2(-vkdiv.*sin(thetamn))); 
    Su = vav(1:mn,1:nn); 
    vav = real(ifft2(vkdiv.*cos(thetamn))); 
    Sv = vav(1:mn,1:nn);
end


