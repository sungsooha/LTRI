function osa = cbct_osaLineRectAngle_v2(angle, dist, rectsz)

if (dist==0)
  osa = 0.5*rectsz(1)*rectsz(2);
  return;
end

if angle==90 || angle==180 
  osa = osaSpecialCase(angle, -dist, rectsz);
elseif angle==0 || angle==270 || angle==360
  osa = osaSpecialCase(angle, dist, rectsz);  
else
  osa = osaGeneralCase_v2(angle, dist, rectsz);
end

end