@ECHO OFF
set prog="../bin/BP.exe"

set root="./fbp_slp/"
set input="slp512_360_filtered.mat"

::set root="./fp_cube_object/offcen/"
::set input="cuboid_offcen.mat"

set LL=0
set LR=1
set LD=2
set TR=4
set TT=3

set alutfn="lut/pos/pos_alut_1501_51.mat"
set hlutfn="lut/pos/pos_hlut_1501_26_8.mat"

::call %prog% -r %root% -i %input% -m %TT%
::call %prog% -r %root% -i %input% -m %TR%
call %prog% -r %root% -i %input% -m %LL% -a %alutfn% -v %hlutfn% -b 2	
::call %prog% -r %root% -i %input% -m %LR% -a %alutfn% -v %hlutfn% -b 2
::call %prog% -r %root% -i %input% -m %LD% -a %alutfn% -v %hlutfn% -b 2



