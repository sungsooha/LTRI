function osa = osaGeneralCase_v2(angle, dist, rectsz)

deg2rad = pi/180;
S = rectsz(1)*rectsz(2);
mxP = rectsz(2);
mnP = -mxP;

a = tan(angle*deg2rad);
b = -dist * sqrt(1+a*a);
  
if angle<=45
  % 1.intercept with y-axis
  b = -b;
  P = b;
  if P<=mnP
    osa = S;
  elseif P>=mxP
    osa = 0;
  else
    osa = S - osaLineRect(a, b, rectsz);
  end
elseif angle<=90
  % 2.intercept with x-axis
  b = -b;
  P = -b/a;
  if P<=mnP
    osa = 0;
  elseif P>=mxP
    osa = S;
  else
    osa = S - osaLineRect(a, b, rectsz);
  end
elseif angle<=135
  % 3.intercept with x-axis
  P = -b/a;
  if P<=mnP
    osa = 0;
  elseif P>=mxP
    osa = S;
  else
    osa = osaLineRect(a, b, rectsz);
  end
elseif angle<=180
  % 4. intercept with y-axis
  P = b;
  if P<=mnP
    osa = 0;
  elseif P>=mxP
    osa = S;
  else
    osa = osaLineRect(a, b, rectsz);
  end
elseif angle <= 225
  % 5. intercept with y-axis
  P =b;
  if P<=mnP
    osa = 0;
  elseif P>=mxP
    osa = S;
  else
    osa = osaLineRect(a, b, rectsz);
  end
elseif angle <= 270
  % 6. intercept with x-axis
  P = -b/a;
  if P<=mnP
    osa = S;
  elseif P>=mxP
    osa = 0;
  else
    osa = osaLineRect(a, b, rectsz);
  end
elseif angle <= 315
  % 7. intercept with x-axis
  b = -b;
  P = -b/a;
  if P<=mnP
    osa = S;
  elseif P>=mxP
    osa = 0;
  else
    osa = S - osaLineRect(a, b, rectsz);
  end
else
  %8. intercept with y-axis
  b = -b;
  P = b;
  if P<=mnP
    osa = S;
  elseif P>=mxP
    osa = 0;
  else
    osa = S - osaLineRect(a, b, rectsz);
  end  
end

end











