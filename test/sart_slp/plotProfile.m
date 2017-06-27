function plotProfile
clear all;
addpath('../PlotPub/lib');
load('slp512_360.mat','nxyz', 'dxyz');
load('slice_im.mat','xy', 'zx', 'yz');
load('slice_ll.mat','ll_xy', 'll_zx', 'll_yz');
load('slice_lr.mat','lr_xy', 'lr_zx', 'lr_yz');
load('slice_ld.mat','ld_xy', 'ld_zx', 'ld_yz');
load('slice_tr.mat','tr_xy', 'tr_zx', 'tr_yz');
load('slice_tt.mat','tt_xy', 'tt_zx', 'tt_yz');

xx = ((1:nxyz(1))-1 - 0.5*(nxyz(1)-1))*dxyz(1);
yy = ((1:nxyz(2))-1 - 0.5*(nxyz(2)-1))*dxyz(2);
zz = ((1:nxyz(3))-1 - 0.5*(nxyz(3)-1))*dxyz(3);


% figure;
% plot(yy, ...
%   [xy(256,:)',...
%    ll_xy(256,:)',lr_xy(256,:)',ld_xy(256,:)',...
%    tr_xy(256,:)',tt_xy(256,:)']);
% opt = plotOption([yy(end),yy(1)],'./img/xy_v.png','North');
% setPlotProp(opt);

% figure;
% iy = 410;
% plot(xx, ...
%   [xy(:,iy),...
%    ll_xy(:,iy),lr_xy(:,iy),ld_xy(:,iy),...
%    tr_xy(:,iy),tt_xy(:,iy)]);
% opt = plotOption([xx(1),xx(end)],'./img/xy_h.png','North');
% setPlotProp(opt);

% figure;
% ix = 193;
% plot(zz, ...
%   [zx(ix,:)',...
%    ll_zx(ix,:)',lr_zx(ix,:)',ld_zx(ix,:)',...
%    tr_zx(ix,:)',tt_zx(ix,:)']);
% opt = plotOption([zz(1),zz(end)],'./img/zx_v.png','North');
% setPlotProp(opt);

% figure;
% iy = 256;
% plot(zz, ...
%   [zx(:,iy),...
%    ll_zx(:,iy),lr_zx(:,iy),ld_zx(:,iy),...
%    tr_zx(:,iy),tt_zx(:,iy)]);
% opt = plotOption([zz(1),zz(end)],'./img/zx_h.png','South');
% setPlotProp(opt);

% figure;
% ix = 283;
% plot(zz, ...
%   [yz(ix,:)',...
%    ll_yz(ix,:)',lr_yz(ix,:)',ld_yz(ix,:)',...
%    tr_yz(ix,:)',tt_yz(ix,:)']);
% opt = plotOption([zz(1),zz(end)],'./img/yz_v.png','North');
% setPlotProp(opt);

figure;
iy = 192;
plot(yy, ...
  [yz(:,iy),...
   ll_yz(:,iy),lr_yz(:,iy),ld_yz(:,iy),...
   tr_yz(:,iy),tt_yz(:,iy)]);
opt = plotOption([yy(end),yy(1)],'./img/yz_h.png','North');
setPlotProp(opt);

end

function opt = plotOption(xlim, fn, lloc)
opt.XLabel = 'distance';
opt.YLabel = 'density (HU)';
opt.BoxDim = [5, 5];
%opt.XLim = [xx(1),xx(end)];
opt.XLim = xlim;
%opt.XLim = [zz(1),zz(end)];
opt.Colors = [
  0, 0, 0;
  1, 0, 0;
  1, 0.25, 1;
  0, 0, 1;
  0, 1, 0;
  0.5, 1, 1.  
];

opt.LineWidth = [2,2];
opt.LineStyle = {'-',':','--','-.',':','--'};
opt.Markers = {'none', 'o', 's', '^', '>', '<'};
opt.MarkerSpacing = [0, 150, 170, 190, 210, 230];
opt.Legend={'SLP','LL','LR','LD','TR','TT'};

opt.LegendBox = 'on';%    bounding box of legend: 'on'/'off'; default: 'off'
opt.LegendBoxColor = [1 1 1]; %: color of the bounding box of legend; default: 'none'
%   LegendTextColor: color of the legend text; default: [0,0,0]
opt.LegendLoc = lloc;%:    'NorthEast', ..., 'SouthWest': legend location

opt.FileName = fn;
%opt.FileName = 'profile_xy_vertical.png';
%opt.FileName = 'profile_yz_hor.png';
end