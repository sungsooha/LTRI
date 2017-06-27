function osa = osaLineSquare(ang, p, side)
% ang: ang <= 0 <= 45, [in deg]
% py: intercept point with y-axis
% side: side length of a square

if p==0, osa = 0.5*side*side; return, end
if ang> 45 && ang<=90 , osa = osaLineSquare(90 - ang, p, side); return, end
if ang> 90 && ang<=135, osa = osaLineSquare(ang - 90, p, side); return, end
if ang>135 && ang<=180, osa = side^2 - osaLineSquare(180 - ang, p, side); return, end
if ang>180 && ang<=225, osa = side^2 - osaLineSquare(ang - 180, p, side); return, end
if ang>225 && ang<=270, osa = side^2 - osaLineSquare(270 - ang, p, side); return, end
if ang>270 && ang<=360, osa = side^2 - osaLineSquare(ang - 270, p, side); return, end

% t = sn*x + cs*y
ang = deg2rad(ang);
cs = cos(ang);
sn = sin(ang);
t = p*cs;

s2 = 0.5*side;
v = [
             -s2,              s2;
  (t - s2*cs)/sn,              s2;
              s2,              s2;
              s2,  (t - s2*sn)/cs;
              s2,             -s2;
  (t + s2*cs)/sn,             -s2;
             -s2,             -s2;
             -s2,  (t + s2*sn)/cs];

check = sum(isinf(v),2) | sum(isnan(v),2) | sum(abs(v)>s2,2);
        
v = v( ~check , : );
v = v( v(:,1)*sn + v(:,2)*cs - t>= -1e-15, :);
v = unique(v,'rows','stable');

osa = polyarea(v(:,1), v(:,2));



end