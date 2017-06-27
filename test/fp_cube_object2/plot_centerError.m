function plot_centerError
addpath('../PlotPub/lib');
load center_error.mat

angles = 0:90;
d_ll = d_ll(1:91)*0.1;
d_lr = d_lr(1:91)*0.1;
d_ld = d_ld(1:91)*0.1;
d_tt = d_tt(1:91)*0.1;
d_tr = d_tr(1:91)*0.1;

figure;
plot(angles, d_ll);
hold on;
plot(angles, d_lr);
plot(angles, d_ld);
plot(angles, d_tt);
plot(angles, d_tr);
hold off

opt.XLabel = 'Projection angles (degree)';
opt.YLabel = 'Max. Abs. difference';
opt.BoxDim = [5, 5];
opt.XLim = [0,90];
opt.Colors = [
  1, 0, 0;
  1, 0.25, 1;
  0, 0, 1;
  0, 1, 0;
  0, 0, 0
];

opt.YScale = 'log';

opt.LineWidth = [2,2,2,2,2];
opt.LineStyle = {'-',':','--','-.','-'};
opt.Markers = {'o', 's', '^', '>', '<'};
opt.MarkerSpacing = [8, 10, 15, 15, 9];
opt.Legend={'LL','LR','LD','TT','TR'};

opt.LegendBox = 'on';%    bounding box of legend: 'on'/'off'; default: 'off'
opt.LegendBoxColor = [1 1 1]; %: color of the bounding box of legend; default: 'none'
%   LegendTextColor: color of the legend text; default: [0,0,0]
opt.LegendLoc = 'SouthEast';%:    'NorthEast', ..., 'SouthWest': legend location

opt.FileName = 'centerError.png';

setPlotProp(opt);
end










