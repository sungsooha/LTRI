function vol = cbct_osvPlaneCuboid(n, d, volsz)

  n = double(n);
  d = double(d);
  volsz = double(volsz);

  vol = 0;

  vertex = 0.5*[-1, -1, -1; -1,  1, -1; 1,  1, -1; 1, -1, -1;
                -1, -1,  1; -1,  1,  1; 1,  1,  1; 1, -1,  1 ];
  vertex = vertex.*repmat(volsz, size(vertex,1), 1);

  dist = dot(vertex, repmat(n, size(vertex,1), 1), 2) + d;
  h = dist>=0;

  if sum(h)==0
    vol = 0;
  elseif sum(h)==8
    vol = volsz(1)*volsz(2)*volsz(3);
  else
    
    edge = [ 1, 2; 2, 3; 3, 4; 4, 1;
             5, 6; 6, 7; 7, 8; 8, 5;
             1, 5; 2, 6; 3, 7; 4, 8 ];
    pset = zeros(12,3);
    np = 0;
    for e=1:12
      [p, flag] = pLinePlane(n,d,vertex(edge(e,1),:),vertex(edge(e,2),:));
      if flag
        np = np+1;
        pset(np,:) = p;
      end
    end
    pset = pset(1:np,:);
    vertex = vertex(h==1,:);
    P = [pset;vertex];
    DT = delaunayTriangulation(P);
    if ~isempty(DT.ConnectivityList)
      vol = volTetrahedron(DT.Points, DT.ConnectivityList);
    end
  end
  
  vol = single(vol);
end