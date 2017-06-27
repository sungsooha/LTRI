@ECHO OFF
set prog="../bin/FP.exe"

::set root="./fp_cube_object2/center/"
::set input="proj_cuboid_center.mat"
::set haveOffset=0

set root="./fp_cube_object2/offcen/"
set input="proj_cuboid_offcen.mat"
set haveOffset=1

set LL=0
set LR=1
set LD=2
set TR=4
set TT=3

set alutfn="../lut/pos/alut_pos_1501_51.mat"
set hlutfn="../lut/pos/hlut_pos_1501_26_8.mat"

call %prog% -r %root% -i %input% -m %TT% -ost %haveOffset%
call %prog% -r %root% -i %input% -m %TR% -ost %haveOffset%
call %prog% -r %root% -i %input% -m %LL% -a %alutfn% -v %hlutfn% -ost %haveOffset%	
call %prog% -r %root% -i %input% -m %LR% -a %alutfn% -v %hlutfn% -ost %haveOffset%	
call %prog% -r %root% -i %input% -m %LD% -a %alutfn% -v %hlutfn% -ost %haveOffset%	



