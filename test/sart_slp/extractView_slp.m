function [xy, zx, yz] = extractView_slp(slp, nxyz, dxyz, clim)

xx = ((1:nxyz(1))-1 - 0.5*(nxyz(1)-1))*dxyz(1);
yy = ((1:nxyz(2))-1 - 0.5*(nxyz(2)-1))*dxyz(2);
zz = ((1:nxyz(3))-1 - 0.5*(nxyz(3)-1))*dxyz(3);

% xy, 
figure,
xy = slp(:,:,190);
im(xx, yy, xy, clim), cbar

%% zx, 
slpx = rotVolX(slp);
zx = slpx(:,:,256);
figure, im(zz,xx,zx,clim,'cbar');
clear slpx

%%
slpy = rotVolY(slp);
yz = slpy(:,:,224);
figure, im(yy,zz,yz,clim,'cbar');
clear slpy

end