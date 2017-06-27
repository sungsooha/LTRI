function osa = osaLineRect(slop, intercept, rectsz)
if intercept==0
  osa = 0.5*rectsz(1)*rectsz(2);
  return;
end

w2 = 0.5*rectsz(1);
h2 = 0.5*rectsz(2);

% find vertices that satisfy slop*x - y + b >= 0
v = npermutek([-1, 1], 2);
v = v.*repmat([w2, h2], size(v,1), 1);
check = dot(v, repmat([slop, -1], size(v,1), 1), 2) + intercept;
v = v(check>0,:);

% find intersection points with the given rectangle
p = [ ( h2 - intercept)/slop,                  h2;
                          w2, slop*w2 + intercept;
      (-h2 - intercept)/slop,                 -h2;
                         -w2, -slop*w2 + intercept ];
check = zeros(4,1);
for i=1:4
  pix = p(i,1);
  piy = p(i,2);
  checkx = -w2<=pix && pix<=w2;
  checky = -h2<=piy && piy<=h2;
  check(i) = checkx && checky;
end
p = p(check>0,:);

vertices = [v;p];
N = size(vertices,1);
if N==1 || N==0
  osa = 0;
else%if (3<=N && N<=5)
  vertices = sortCCW(vertices);
  osa = polyarea(vertices(:,1), vertices(:,2));
% else
%   error('cbct_osaLineRect:vertex',...
%     'unexpected number of vertices is found %d', N);
end
end