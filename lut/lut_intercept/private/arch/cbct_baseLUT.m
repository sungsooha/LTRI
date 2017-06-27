function bLUT = cbct_baseLUT(angle, intercept, side)
rectsz = [side, side];
nangle = length(angle);
nintercept = length(intercept);

data = zeros(nangle, nintercept);
total = nangle*nintercept;
count = 0;

hProgress = waitbar(0, 'building bLUT: 0%');
for j=1:nangle
  for i=1:nintercept
    data(j,i) = cbct_osaLineRectAngle(angle(j), intercept(i), rectsz);

    count = count + 1;
    perc = count/total;
    waitbar(perc,hProgress,sprintf('building bLUT: %.1f%%',perc*100));
  end
end
close(hProgress);

bLUT = struct(...
  'X', angle, 'Y', intercept, 'V', data,...
  'side', side, 'S', side*side);
end






