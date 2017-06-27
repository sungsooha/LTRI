function [xy, zx, yz] = extractView(xyz, nxyz, dxyz, clim)

xx = ((1:nxyz(1))-1 - 0.5*(nxyz(1)-1))*dxyz(1);
yy = ((1:nxyz(2))-1 - 0.5*(nxyz(2)-1))*dxyz(2);
zz = ((1:nxyz(3))-1 - 0.5*(nxyz(3)-1))*dxyz(3);

% xy, 
figure,
xy = xyz(:,:,193); % x, y, z
im(xx, yy, xy, clim), cbar

%% zx, 
%slpx = rotVolX(vol);
zxy = permute(xyz,[3 1 2]); % y z x
zx = zxy(:,:,222);
figure, im(zz,xx,zx,clim,'cbar');
%clear slpx

%% yz
yzx = permute(zxy, [3 1 2]);
%slpy = rotVolY(vol);
yz = yzx(:,:,247);
figure, im(yy,zz,yz,clim,'cbar');
%clear slpy

end