% Batch evaluate

%% Clean workspace
clear;
clc;
close all;

%% Add 'process' path

addpath('../Process_scripts/');
addpath('../Analysis/');

%% Load forcing details

% forceName = '_input_forcing_SMgeo_5day.mat';
% % forceName = '_input_forcing_2DQGn.mat';
% forcePath = 'forcings/';
% if ~exist('forcing','var')
%     forcing = readForcing(forceName,forcePath);
% end

%% Load data

batchCases = {'FSD_300_7500_2_geo','FSD_300_2500_2_geo','FSD_300_1000_2_geo'};%'FSD_3000_30000_2_geo'};%,'FSD_300_1000_2_geo'};
% batchCases = {'FSD_3000_30000_2_geo'};

B = length(batchCases);
bins = logspace( log10(1e2), log10(10e3), 10 );
% bins = logspace( log10(3000), log10(30e3), 12 );
% bins = [ logspace( log10(300), log10(30e3), 20 ) ];
logBinSz = 10.^binCenters(log10(bins));

cols = lines(B);
% Hs = ["f49097","e0ca3c","577590","2ebfa5","5f0f40"];
% Hs = ["dc0073","ffd166","1b998b","2191fb","ff5714"];
% cols = hex2colour(Hs);

for b = 1:B

    fName = strcat('floeOut_',batchCases{b} );
    load( strcat('../MAT_Files/',fName) );
    

    
    N = length(floeOut);
    SIC = arrayfun( @(S) S.chunk.sic, floeOut ); 
    idL = arrayfun( @(S) mean(S.chunk.ioblDiss(2:end)), floeOut );
    idSL = arrayfun( @(S) mean(S.chunk.ioblSteadyDiss(2:end)), floeOut );

    intDiss = NaN(19,N);


    binDiss = NaN( length(bins)-1,N );  
    binSDiss = NaN( length(bins)-1,N );  
    binEDiss = NaN( length(bins)-1,N );  
    for n = 1:N
 
        % Extract data
        [M,L] = size(floeOut(n).state.ui);
        SIC(n) = floeOut(n).chunk.sic;
        ind = 2:L;        
        dissML = floeOut(n).ocean.avg.diss(:,ind);
        sDissML = floeOut(n).ocean.avg.steadyDiss(:,ind);
        eDissML = floeOut(n).ocean.avg.excessDiss(:,ind);
%         eDissML = dissML - sDissML;
        
       comUo = floeOut(n).ocean.com.uo(:,ind)+1i*floeOut(n).ocean.com.vo(:,ind);
%         gInd = ( abs(comUo)>0.02 );
%         gInd = ( sDissML>1e-6 );
        gInd = 1:numel(dissML);
        floeSizeMat = floeOut(n).floes.floeSize.'*ones(1,size(dissML,2));
    
        % Bin by size
        % (median values)
        avgFun = @mean;
        [binSz,binDiss(:,n)]  = binFun( floeSizeMat(gInd), dissML(gInd),  bins,binFunction=avgFun );
        [binSz,binEDiss(:,n)] = binFun( floeSizeMat(gInd), eDissML(gInd), bins,binFunction=avgFun );
        [binSz,binSDiss(:,n)] = binFun( floeSizeMat(gInd), sDissML(gInd), bins,binFunction=avgFun );
        
        intDiss(2:L,n) = floeOut(n).chunk.ioblDiss(2:L);
    end
    
   %% Plot bin-averaged dissipation

    fH = figure(3); 
    fH.Units = 'inches';
    fH.Position([3,4]) = [12,3];
    fontSize = 12;

%     ax = gca;
    ax = subplot(1,B,b);
    hold on;
%     pH = plot( binSz, binSDiss,'s-',MarkerSize=6, LineWidth=1.5);
    pH = plot( logBinSz, binEDiss,'d-',MarkerSize=6, LineWidth=1.5);
%     pH = plot( binSz, binDiss,'o-',MarkerSize=6, LineWidth=1.5,MarkerFaceColor='w');
    cShift = (0.9)*(2*SIC-1);
    lineCols =  colourShift(cols(b,:),cShift);
    [pH.Color] = disperse( lineCols.' );


    ax.XLim = [300,1e4];    ax.XTick = 10.^(2:5);
    ax.YLim = [8e-6,2e-3];  ax.YTick = 10.^(-6:0);
%     ax.YLim = [1e-8,2e-3];
%     ax.YLim = [1e0,2e5];
    
    ax.XScale = 'log';
    ax.YScale = 'log';
    ax.FontSize = fontSize;
    grid on;

    ylabel({'Mean floe-averaged';'excess disspation [W/m^2]'});
    xlabel('Floe size [m]');

%     legLab = sprintfc('%2.0f%%',100*SIC);
%     legend(pH,legLab,NumColumns=2,Location='southeast');


    
    %% Plot domain-integrated dissipation

    fH = figure(4);
    fH.Units = 'inches';
    fH.Position([3,4]) = [6,3];
      
    ax = gca;    
    hold on;
    pH5(b) = plot( SIC, idL, 'o-', LineWidth=1.5, Color=cols(b,:),MarkerFaceColor='w');
    plot( SIC, idSL, '--', LineWidth=1, Color=cols(b,:) );

    ax.XLim = [0,1];
    ax.YLim = [5e3,1e7];
    ax.YScale = 'log';
    ax.FontSize = fontSize;
    grid on;

    ylabel({'Integrated disspation [W]'});
    xlabel('Sea ice concentration');

%     if b==B
%     legend( pH5 , {'big floes','medium floes','small floes'},...
%            Location='southeast')
%     end

    %% Plot domain-integrated dissipation timeseries

    fH = figure(5); 
    fH.Units = 'inches';
    fH.Position([3,4]) = [12,3];
    fontSize = 12;

    T = (0:18)*4*3600/86400;
%     ax = gca;
    ax = subplot(1,B,b);
    hold on;
    pH = plot( T, intDiss, LineWidth=1.5);
    cShift = (0.9)*(2*SIC-1);
    lineCols =  colourShift(cols(b,:),cShift);
    [pH.Color] = disperse( lineCols.' );
    ax.YLim = [5e3,1e7];
    ax.YScale = 'log';
    ax.FontSize = fontSize;
    ax.XLim = [0,3];
    grid on;

    xlabel('Time [days]');
    ylabel('Integrated dissipation [W]')

end




%%
%{
fH = figure(10); clf;
fH.Units = 'inches';
fH.Position([3,4]) = [6,3];

ax = gca;
hold on;

edispos = eDissML>0;
scatter( floeSizeMat(edispos), eDissML(edispos),9,cols(b,:),'filled',MarkerFaceAlpha=0.85 );
% scatter( floeSizeMat(~edispos),-eDissML(~edispos),9,'r','filled',MarkerFaceAlpha=0.85 );
% plot( binSz, binEDiss,'ko-',MarkerSize=6, LineWidth=1.5);

ax.XLim = [300,1e4];
% ax.YLim = [1e-10,1e0];
ax.XScale = 'log'; 
ax.YScale = 'log';
ax.FontSize = fontSize;

grid on;
xlabel('Floe size [m]');
ylabel({'floe-averaged steay disspation';'[W/m^2]'});
%}