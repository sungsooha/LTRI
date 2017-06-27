function osv = osvPlaneCuboid(n, d, voxsz)

if d==0, osv = 0.5*voxsz(1)*voxsz(2)*voxsz(3); return, end

osv = 0;
vertex = 0.5*[-1, -1, -1; -1,  1, -1; 1,  1, -1; 1, -1, -1;
              -1, -1,  1; -1,  1,  1; 1,  1,  1; 1, -1,  1 ];
vertex = vertex.*repmat(voxsz, size(vertex,1), 1);

dist = dot(vertex, repmat(n, size(vertex,1), 1), 2) - d;
h = dist>=0;

if sum(h)==0, osv = 0; return, end
if sum(h)==8, osv = voxsz(1)*voxsz(2)*voxsz(3); return, end
    
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
P = unique(P,'rows','stable');

DT = delaunayTriangulation(P);
if ~isempty(DT.ConnectivityList)
  osv = volTetrahedron(DT.Points, DT.ConnectivityList);
end

end