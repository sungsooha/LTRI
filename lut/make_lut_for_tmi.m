%load('ph_2.mat','dxyz', 'dsd', 'nst', 'dst');

% dxyz = [1, 1, 1];
% dsd =  949.075;
% nst = [888, 640];
% dst = [1.0279, 1.0964];

% dxyz = [0.415, -0.415, 0.83];
% dsd =  1147;
% nst = [1024, 384];
% dst = [0.3842, 0.7558];

bBuild_hlut = true;
bBuild_alut = false;

if (bBuild_hlut)
  dxyz = abs(dxyz);
  dx = dxyz(1);
  dy = dxyz(2);
  dz = dxyz(3);
  nd = 4000;
  na = 45;
  np = 60;
  maxdist = 0.5*sqrt(dot(dxyz,dxyz));
  maxazi = 45;
  maxpol = ceil(rad2deg(atan( (0.5*nst(2)*dst(2))/dsd ))) + 5;

  [lut, ddist, dazi, dpol] = ...
    build_lut_height(nd, na, np, maxdist, maxazi, maxpol, dx, dy, dz);

  [nd, np, na] = size(lut);
  %lut = permute(lut, [2 1 3]);
  lut = single(lut);
  
  
  save(['hlut_' num2str(nd) '_' num2str(na) '_' num2str(np) '.mat'],...
    'lut', 'nd', 'np', 'na', 'ddist', 'dazi', 'dpol', 'dx', 'dy', 'dz');
end

if (bBuild_alut)
  dx = dxyz(1);
  dy = abs(dxyz(2));
  nd = 4000;
  na = 90;
  maxdist = 0.5*sqrt(dx*dx + dy*dy);
  maxang = 45;
  [lut, ddist, dang] = build_lut_area(nd, na, maxdist, maxang, dx, dy);

  [nd, na] = size(lut);

  %lut = permute(lut, [2 1]);
  lut = single(lut);
  save(['alut_' num2str(nd) '_' num2str(na) '.mat'], ...
    'lut', 'nd', 'na', 'ddist', 'dang', 'dx', 'dy');
end

clear all
