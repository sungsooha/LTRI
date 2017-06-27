function osa = cbct_osaLineRectAngle(angle, intercept, rectsz)
if ~((0<=angle) && (angle<=360))
  error('cbct_osaLineRectAngle:input', '0<=angle<=360');
end

if (intercept==0)
  osa = 0.5*rectsz(1)*rectsz(2);
  return;
end

if angle==0 || angle==90 || angle==180 || angle==270
  osa = osaSpecialCase(angle, intercept, rectsz);
else
  osa = osaGeneralCase(angle, intercept, rectsz);
end
end