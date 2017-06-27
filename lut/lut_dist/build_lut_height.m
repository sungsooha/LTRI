function [lut, ddist, dazi, dpol] = ...
  build_lut_height(nd, na, np, maxdist, maxazi, maxpol, dx, dy, dz)

if (nargin == 1 && strcmp(nd,'test')), test_lut_height, return, end

volsz = [dx, dy, dz];
base = dx*dy;

ddist = maxdist / nd;
dazi = maxazi / na;
dpol = maxpol / np;

distSet = (0:nd)*ddist;
aziSet = (0:na)*dazi;
polSet = (0:np)*dpol;

lut = zeros(nd+1, np+1, na+1);

mxCount = (na+1)*(np+1);
count = 1;

h = waitbar(0,'Please wait...');
for ia=1:length(aziSet)
  azi = deg2rad(aziSet(ia));
  for ip=1:length(polSet)
    pol = deg2rad(polSet(ip));
    n = [sin(pol)*sin(azi), -sin(pol)*cos(azi), cos(pol)];
    for id=1:length(distSet)
      dist = distSet(id);
      vol = cbct_osvPlaneCuboid(n, -dist, volsz);
      lut(id,ip,ia) = vol / base;
    end
    count = count + 1;
    waitbar(count / mxCount)
  end
  
end
close(h);

end

function test_lut_height

dx = 1; dy = 1; dz = 1;
nd = 1500; na = 2; np = 2;
maxdist = 0.5*sqrt(dx*dx + dy*dy + dz*dz);
maxazi = 45;
maxpol = 45;

[lut, ddist, dazi, dpol] = ...
  build_lut_height(nd, na, np, maxdist, maxazi, maxpol, dx, dy, dz);

distSet = (0:nd)*ddist;


figure,
subplot(131), plot(distSet, [lut(:,1,1),lut(:,2,1),lut(:,3,1)]); % azi = 0, pol = 0
subplot(132), plot(distSet, [lut(:,1,2),lut(:,2,2),lut(:,3,2)]); % azi = 0, pol = 0
subplot(133), plot(distSet, [lut(:,1,3),lut(:,2,3),lut(:,3,3)]); % azi = 0, pol = 0


end