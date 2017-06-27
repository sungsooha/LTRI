function build_hlut_pos(npos, npol, nazi, voxsz,...
  maxAzi, maxPos, maxPol, fn)

%maxAzi = 45;
%maxPos = voxsz(3);
%maxPol = 45;

step_azi = maxAzi/nazi;
step_pos = maxPos/npos;
step_pol = maxPol/npol;

aziSet = (0:nazi)*step_azi;
posSet = (0:npos)*step_pos;
polSet = (0:npol)*step_pol;

npos = length(posSet);
nazi = length(aziSet);
npol = length(polSet);

hlut = zeros(npos,npol,nazi);

h = waitbar(0,'Please wait...');
mxCount = nazi*npol;
count = 0;

for ia=1:nazi
  azi = deg2rad(aziSet(ia));
  for ip=1:npol
    pol = deg2rad(polSet(ip));
    n = [sin(pol)*sin(azi), -sin(pol)*cos(azi), cos(pol)];
    for pp=1:npos
      d = n(3)*posSet(pp);
      vol = osvPlaneCuboid(n,d,voxsz);
      hlut(pp,ip,ia) = vol / (voxsz(1)*voxsz(2));
    end
    count = count + 1;
    waitbar(count / mxCount)    
    
  end  
end
close(h);

hlut = single(hlut);

if nargin==8
  save([fn '_' num2str(npos) '_' num2str(npol) '_' num2str(nazi) '.mat'],...
    'hlut', 'npos', 'npol', 'nazi', 'maxAzi', 'maxPos', 'maxPol', ...
    'voxsz', 'step_azi', 'step_pos', 'step_pol');
end

end