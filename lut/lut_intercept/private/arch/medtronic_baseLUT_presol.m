function medtronic_baseLUT_presol
%% base look-up table
scale = 1;
side = 0.415;

presol = 400;
viewresol = 360;

step_view = 360/viewresol;
step_pos = (2*side)/presol;

viewangle = step_view*(0:viewresol);
intercept = -(step_pos*(0:presol) - side);

bLUT = cbct_baseLUT(viewangle, intercept, side);

V = bLUT.V;
%V = fliplr(V);

filename = 'baseLUT_360_p400_a360.txt';
fid = fopen(filename,'w');
fprintf(fid, '%d %d\n', length(intercept), length(viewangle));
fprintf(fid, '%f %f %f\n', step_pos, step_view, scale);
for j=1:length(viewangle)
  for i=1:length(intercept)
    fprintf(fid,'%f', scale*V(j,i));
    if i<length(intercept)
      fprintf(fid,' ');
    end
  end
  if j<length(viewangle)
    fprintf(fid,'\n');
  end
end
fclose(fid);

dlmwrite('intercept_360_p400_a360.txt', intercept, ' ');

end