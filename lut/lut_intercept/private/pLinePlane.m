function [p, flag] = pLinePlane(n,d,p0,p1)
  p = 0;
  flag = false;

  d0 = dot(p0,n) - d;
  d1 = dot(p1,n) - d;

  if d0*d1<0
    d = d0 / (d0-d1);
    p = p0 + d*(p1-p0);
    flag = true;
  end
end