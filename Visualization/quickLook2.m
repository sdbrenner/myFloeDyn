
clear;
clc;
% close all;


outName = 'out_Sp7jZ.h5';
fo = readFloeOut(outName);
[M,L] = size(fo.state.xa);
foL = subsetFloeOut(fo,L);
foL = makeDomainPoly(foL);

fH = figure(1); %clf; 
fH.Units = 'inches';
fH.Position([3,4]) = [12,6];

ax = subplot(1,2,1); cla;
plot( fo.state.xa.', fo.state.ya.' );
ax.XLim = fo.chunk.win(1:2);
ax.YLim = fo.chunk.win(3:4);
daspect([1,1,1]);


ax = subplot(1,2,2);
hold on;
plot( foL.chunk.icePoly{end} );
ax.XLim = fo.chunk.win(1:2);
ax.YLim = fo.chunk.win(3:4);
daspect([1,1,1]);



