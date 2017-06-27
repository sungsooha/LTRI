function sorted = sortCCW(unsorted)
N = size(unsorted,1);
angles = zeros(N,1);

for i=1:N
  vec = unsorted(i,:);
  nvec = vec/norm(vec);
  rad = acos(nvec(1));
  if (nvec(2)>=0)
    angles(i) = rad*180/pi;
  else
    angles(i) = 360 - rad*180/pi;
  end
end

[~,index] = sort(angles,'ascend');

sorted = zeros(N,2);
for i=1:N
  sorted(i,:) = unsorted(index(i),:);
end
end