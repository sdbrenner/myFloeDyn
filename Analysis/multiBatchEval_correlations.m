% Batch evaluate

%% Clean workspace
clear;
clc;
close all;

%% Add 'process' path

addpath('../Process_scripts/');
addpath('../Analysis/');

%% Load data

batchCases = {'FSD_300_7500_2_geo','FSD_300_2500_2_geo','FSD_300_1000_2_geo'};

for b = 1:length(batchCases)

    fName = strcat('floeOut_',batchCases{b} );
    load( strcat('../MAT_Files/',fName) );
    
    %% Calculate correlations
    
    N = length(floeOut);
    R = NaN( 2,N );
    SIC = NaN( 1,N );
    
    for n = 1:N
        SIC(n) = floeOut(n).chunk.sic;
        Ui = floeOut(n).state.ui  + 1i*floeOut(n).state.vi;
        Uo = floeOut(n).ocean.avg.uo + 1i*floeOut(n).ocean.avg.vo;
    %     Uiss = seaIceSteadyState( (5e3/(1.4e-4*1)), Uo );
    
    %     rr = corrcoef( abs(Uo(:,2:end)), abs(Ui(:,2:end)) );
        rr = corrcoef( (Uo(:,2:end)), (Ui(:,2:end)) );
        R(1,n) = abs(rr(1,2));
    %     rr = corrcoef( real(Uiss(:,2:end)), real(Ui(:,2:end)) );
    %     R(2,n) = rr(1,2);
    %     rr = corrcoef( imag(Uiss(:,2:end)), imag(Ui(:,2:end)) );
    %     R(3,n) = rr(1,2);
    
        rr = corrcoef( floeOut(n).state.om(:,2:end), floeOut(n).ocean.avg.zeta(:,2:end) );
        R(4,n) = rr(1,2);
    end
    
    %% Plot correlations
    
    % fontName = 'Century Gothic';
    fontName = 'Montserrat';
    fontSize = 14;
    mea = 1;
    mfa = 0.5; 
    ms = 12;
    cols = lines( length(batchCases) );
    
    fH = figure(1); 
    fH.Units = 'inches';
    fH.Position([3,4]) = [5.5,4.5];
    hold on;
    colororder( cols([1,1,2,2,3,3],:) );
    
    % ax = subplot(2,1,1);
    ax = gca;
    hold on;
    pH1(b) = plot( SIC, R(1,:), '.-',MarkerSize=ms,LineWidth=2);
    % scatter( SIC, R(1,:),ms^2,'filled',...
    %          MarkerEdgeColor=grey(0.5),MarkerEdgeAlpha=mea, MarkerFaceAlpha=mfa );
    % scatter( SIC, R(2,:),(0.65*ms)^2,'filled',...
    %          MarkerEdgeColor=grey(0.5),MarkerEdgeAlpha=mea, MarkerFaceAlpha=mfa );
    % scatter( SIC, R(3,:),(0.65*ms)^2,'filled',...
    %          MarkerEdgeColor=grey(0.5),MarkerEdgeAlpha=mea, MarkerFaceAlpha=mfa );
    ax.XLim = [0,1];
    ax.YLim = [0,1];
    ax.YTick = 0:0.25:1;
    ax.FontName = fontName;
    ax.FontSize = fontSize;
    grid on;
    xlabel('Sea ice concentration')
    % ylabel('Velocity correlation')
    % ylabel('Corr[ {\bfU}_i, \langle{\bfu}_o\rangle ]')
    
    % ax = subplot(2,1,2);
    hold on;
    pH2(b) = plot( SIC, R(4,:), 'o--',MarkerSize=ms/3,LineWidth=2,MarkerFaceColor='w');
    % scatter( SIC, R(4,:),ms^2,'filled',...
    %          MarkerEdgeColor=grey(0.5),MarkerEdgeAlpha=mea, MarkerFaceAlpha=mfa );
    ax.XLim = [0,1];
    ax.YLim = [0,1];
    ax.YTick = 0:0.25:1;
    ax.FontName = fontName;
    ax.FontSize = fontSize;
    grid on;
    xlabel('Sea ice concentration')
    % ylabel('Rotation correlation')
    % ylabel('Corr[ \omega_i, \langle\zeta_o\rangle ]')
    
    ylabel('Correlation coefficient')


end
%% Create correlations legend

fsdLabels = [ repmat({''},1,length(batchCases)), strrep(batchCases,'_','-') ];%...
%            'FSD: 300—7500, \alpha=2'; 'FSD: 300—2500, \alpha=2'; 'FSD: 300—1000, \alpha=2'};

legFontSize = 10;
legH1 = legend( ax, [pH1;pH2].',fsdLabels);
legH1.Position = [ax.Position(1),ax.Position(2),0.6,0.15];
legH1.FontSize = legFontSize;
legH1.Box = 'off';
legH1.NumColumns = 2;

tH1 = text(0,0,'{\itR}[ {\bfU}_i,\langle{\bfu}_o\rangle ]',...
          Units='normalized',FontName=fontName,FontSize=legFontSize);
tH1.Position = [0.03,0.2];

tH2 = text(0,0,'{\itR}[ \omega_i,\langle\zeta_o\rangle ]',...
          Units='normalized',FontName=fontName,FontSize=legFontSize);
tH2.Position = [0.23,0.2];
