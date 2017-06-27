@ECHO OFF
set prog="../bin/FP.exe"

set root="./fp_cube_object/center/"
set input="cuboid_center.mat"

::set root="./fp_cube_object/offcen/"
::set input="cuboid_offcen.mat"

set LL=0
set LR=1
set LD=2
set TR=4
set TT=3

set alutfn="../lut/pos/pos_alut_8001_181.mat"
set hlutfn="../lut/pos/pos_hlut_2001_21_11.mat"

call %prog% -r %root% -i %input% -m %TT%
call %prog% -r %root% -i %input% -m %TR%
call %prog% -r %root% -i %input% -m %LL% -a %alutfn% -v %hlutfn%	
call %prog% -r %root% -i %input% -m %LR% -a %alutfn% -v %hlutfn%	
call %prog% -r %root% -i %input% -m %LD% -a %alutfn% -v %hlutfn%	



