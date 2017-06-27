function [lut, ddist, dang] = ...
  build_lut_area(nd, na, maxdist, maxang, dx, dy)

if (nargin == 1 && strcmp(nd,'test')), test_lut_area, return, end

ddist = maxdist / nd;
dang = maxang / na;
rectsz = [dx, dy];

angSet = dang*(0:na);
distSet = ddist*(0:nd);

lut = zeros(nd+1, na+1);

h = waitbar(0,'Please wait...');
for j=1:length(angSet)
  ang = angSet(j);
  for i=1:length(distSet)
    dist = distSet(i);
    area = cbct_osaLineRectAngle_v2(ang, dist, rectsz);
    lut(i,j) = area;
  end
  waitbar(j / (na+1))
end
close(h);

end

function test_lut_area

dx = 1; dy = 1;
nd = 1500; na = 50;
maxdist = 0.5*sqrt(dx*dx + dy*dy);
maxang = 45;

[lut, ddist, dang] = build_lut_area(nd, na, maxdist, maxang, dx, dy);

end