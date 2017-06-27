function plotConv
clear all;
addpath('../PlotPub/lib');
load conv.mat

figure;
plot(iter, me_ll);
hold on;
plot(iter, me_lr);
plot(iter, me_ld);
plot(iter, me_tt);
plot(iter, me_tr);
hold off

opt.XLabel = '# iterations';
opt.YLabel = 'RMSE (HU)';
opt.BoxDim = [5, 5];
opt.XLim = [iter(1),iter(end)];
opt.Colors = [
  1, 0, 0;
  1, 0.25, 1;
  0, 0, 1;
  0, 1, 0;
  0, 0, 0
];

%opt.YScale = 'log';

opt.LineWidth = [2,2,2,2,2];
opt.LineStyle = {'-',':','--','-.','--'};
opt.Markers = {'o', 's', '^', '>', '<'};
opt.MarkerSpacing = [100, 150, 200, 250, 300];
opt.Legend={'LL','LR','LD','TT','TR'};

opt.LegendBox = 'on';%    bounding box of legend: 'on'/'off'; default: 'off'
opt.LegendBoxColor = [1 1 1]; %: color of the bounding box of legend; default: 'none'
%   LegendTextColor: color of the legend text; default: [0,0,0]
opt.LegendLoc = 'NorthEast';%:    'NorthEast', ..., 'SouthWest': legend location

opt.FileName = 'plotConv.png';

setPlotProp(opt);










