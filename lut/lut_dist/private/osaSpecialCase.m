function osa = osaSpecialCase(angle, P, rectsz)

rw2 = 0.5*rectsz(1);
rh2 = 0.5*rectsz(2);
base = rectsz(1)*rectsz(2);

if     angle==0
  % P is intercept point with y-axis
  if P<-rh2
    osa = base; 
  elseif P>rh2
    osa = 0;
  else
    osa = (rh2 - P)*rectsz(1);
  end
elseif angle==90
  % P is intercept point with x-axis
  if P<-rw2
    osa = 0; 
  elseif P>rw2
    osa = base;
  else
    osa = (rw2 + P)*rectsz(2);
  end
elseif angle==180
  % P is intercept point with y-axis
  if P<-rh2
    osa = 0; 
  elseif P>rh2
    osa = base;
  else
    osa = (rh2 + P)*rectsz(1);
  end
elseif angle==270
  % P is intercept point with x-axis
  if P<-rw2
    osa = base; 
  elseif P>rw2
    osa = 0;
  else
    osa = (rw2 - P)*rectsz(2);
  end
else  %angle==360 (equal to angle==0)
  % P is intercept point with y-axis
  if P<-rh2
    osa = base; 
  elseif P>rh2
    osa = 0;
  else
    osa = (rh2 - P)*rectsz(1);
  end
end

end