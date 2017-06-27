function build_alut_pos(np, na, side, maxAng, maxPos, fn)

%if nargin == 1 && strcmp(na,'test'), test_build_alut, return, end

%side = 2;
%na = 45;
%np = 1000;

% maxAng = 45;
% maxPos = side;

step_ang = maxAng/na;
step_pos = maxPos/np;

angSet = (0:na)*step_ang;
posSet = (0:np)*step_pos;

na = length(angSet);
np = length(posSet);

alut = zeros(np, na);

h = waitbar(0, 'Please wait...');
for ia = 1:na;
  for ip = 1:np;
    alut(ip,ia) = osaLineSquare(angSet(ia), posSet(ip), side);
  end
  waitbar(ia / na);
end
close(h);

%if nargin==4
  alut = single(alut);
  save([fn '_' num2str(np) '_' num2str(na) '.mat'],...
    'alut', 'na', 'np', 'step_ang', 'step_pos', 'maxAng', 'maxPos', 'side');
%end


end