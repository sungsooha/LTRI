function medtronic_heightLUT_v2
%% height look-up table
scale = 1;

side = 0.415;
height = 0.83;

viewresol = 45;
polarresol = 20;
posresol = 400;

step_view  = 45/viewresol;
step_polar = 10/polarresol;%45/polarresol;
step_pos   = (2*height)/posresol;

viewangles = step_view*(0:viewresol);
polarangle = step_polar*(0:polarresol);
intercept  = step_pos*(0:posresol) - height;

hLUT = cbct_heightLUT(intercept, polarangle, viewangles, [side, side, height]);

fid = fopen('medtronic_heightLUT.txt','w');
fprintf(fid,'%d %d %d\n', length(intercept), length(polarangle), length(viewangles));
fprintf(fid,'%f %f %f %f\n', step_pos, step_polar, step_view, scale);
for k=1:length(viewangles)
  for j=1:length(polarangle)
    for i=1:length(intercept)
      fprintf(fid, '%f', scale*hLUT.V(j,i,k));
      if i<length(intercept)
        fprintf(fid,' ');
      end
    end
    fprintf(fid,'\n');
  end
end

fclose(fid);

end

