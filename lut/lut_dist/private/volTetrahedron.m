function vol = volTetrahedron(vertices, list)
  vol = 0;
  for i=1:size(list,1)
      l = list(i,:);
      a = vertices(l(1),:);
      b = vertices(l(2),:);
      c = vertices(l(3),:);
      d = vertices(l(4),:);

      ad = a - d;
      bd = b - d;
      cd = c - d;

      v = abs(dot(ad, cross(bd,cd))) / 6;
      vol = vol + v;
  end
end