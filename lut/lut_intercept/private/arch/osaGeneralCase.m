function osa = osaGeneralCase(angle, P, rectsz)

deg2rad = pi/180;
S   = rectsz(1)*rectsz(2);
mxP = rectsz(2);
mnP = -mxP;

if         0<angle && angle<=45
  % 1. P is intercept point with y-axis
  if P<=mnP
    osa = S;
  elseif P>=mxP
    osa = 0;
  else
    a = tan(angle*deg2rad);
    b = P;
    osa = S - osaLineRect(a, b, rectsz);
  end
elseif    45<angle && angle<=90
  % 2. P is intercept point with x-axis
  if P<=mnP
    osa = 0;
  elseif P>=mxP
    osa = S;
  else
    a = tan(angle*deg2rad);
    b = -a*P;
    osa = S - osaLineRect(a, b, rectsz);
  end
elseif    90<angle && angle<=135
  % 3. P is intercept point with x-axis
  if P<=mnP
    osa = 0;
  elseif P>=mxP
    osa = S;
  else
    a = tan(angle*deg2rad);
    b = -a*P;
    osa = osaLineRect(a, b, rectsz);
  end
elseif   135<angle && angle<=180
  % 4. P is intercept point with y-axis
  if P<=mnP
    osa = 0;
  elseif P>=mxP
    osa = S;
  else
    a = tan(angle*deg2rad);
    b = P;
    osa = osaLineRect(a, b, rectsz);
  end
elseif   180<angle && angle<=225
  % 5. P is intercept point with y-axis
  if P<=mnP
    osa = 0;
  elseif P>=mxP
    osa = S;
  else
    a = tan(angle*deg2rad);
    b = P;
    osa = osaLineRect(a, b, rectsz);
  end
elseif   225<angle && angle<=270
  % 6. P is intercept point with x-axis
  if P<=mnP
    osa = S;
  elseif P>=mxP
    osa = 0;
  else
    a = tan(angle*deg2rad);
    b = -a*P;
    osa = osaLineRect(a, b, rectsz);
  end
elseif   270<angle && angle<=315
  % 7. P is intercept point with x-axis
  if P<=mnP
    osa = S;
  elseif P>=mxP
    osa = 0;
  else
    a = tan(angle*deg2rad);
    b = -a*P;
    osa = S - osaLineRect(a, b, rectsz);
  end
else   % 315<angle && angle<=360
  % 8. P is intercept point with y-axis
  if P<=mnP
    osa = S;
  elseif P>=mxP
    osa = 0;
  else
    a = tan(angle*deg2rad);
    b = P;
    osa = S - osaLineRect(a, b, rectsz);
  end
end
end