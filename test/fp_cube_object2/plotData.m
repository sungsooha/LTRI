function plotData


%center
root = './center/';
ya = load_proj([root 'proj_cuboid_center.mat']);
fn = 'center_error.mat';
usemax = true;

%offcen
% root = './offcen/';
% ya = load_proj([root 'proj_cuboid_offcen.mat']);
% fn = 'offcen_error.mat';
% usemax = true;

yd_ll = load_proj([root 'fp_LL.mat']);
yd_lr = load_proj([root 'fp_LR.mat']);
yd_ld = load_proj([root 'fp_LD.mat']);
yd_tt = load_proj([root 'fp_TT.mat']);
yd_tr = load_proj([root 'fp_TR.mat']);

if usemax
  d_ll = proj_abs_max(ya, yd_ll);
  d_lr = proj_abs_max(ya, yd_lr);
  d_ld = proj_abs_max(ya, yd_ld);
  d_tt = proj_abs_max(ya, yd_tt);
 d_tr = proj_abs_max(ya, yd_tr);
else
  d_ll = proj_abs_avg(ya, yd_ll);
  d_lr = proj_abs_avg(ya, yd_lr);
  d_ld = proj_abs_avg(ya, yd_ld);
  d_tt = proj_abs_avg(ya, yd_tt);
  d_tr = proj_abs_avg(ya, yd_tr);  
end

%clear ya yd_ll yd_lr yd_ld yd_tt yd_tr root

figure, plot([d_ll,d_lr,d_ld,d_tt,d_tr])

save(fn);
%save center_abs_avg.mat
end

function y = load_proj(matfn)

load(matfn, 'proj');
y = proj;

end

function me = proj_abs_max(ya, yd)

diff = ya - yd;
nslice = size(ya,3);

me = zeros(nslice,1);
for i=1:nslice
  dd = abs(diff(:,:,i));
  %dd(dd>0.001) = 0;
  me(i) = max(abs(dd(:)));
end

end

function me = proj_abs_avg(ya, yd)

diff = ya - yd;
nslice = size(ya,3);

me = zeros(nslice,1);
for i=1:nslice
  dd = diff(:,:,i);
  me(i) = sum( abs(dd(:)) )/nnz(ya(:,:,i));
end

end