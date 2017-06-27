@ECHO OFF
set prog="../bin/sart3d.exe"

set root="./sart_slp/"
set input="slp512_360.mat"

set maxiter=500
set lambda=0.0025
set interval=50
set stop=50


set LL=0
set LR=1
set LD=2
set TR=4
set TT=3

::set alutfn="lut/pos/pos_alut_1501_51.mat"
set alutfn="lut/pos/pos_alut_10001_181.mat"
set hlutfn="lut/pos/pos_hlut_1501_26_8.mat"

::call %prog% -r %root% -i %input% -m %LD% -a %alutfn% -v %hlutfn% -T %interval% -I %maxiter% -L %lambda% -S %stop% -b 2

::call %prog% -r %root% -i %input% -m %TR% -T %interval% -I %maxiter% -L %lambda% -S %stop% -b 2
::call %prog% -r %root% -i %input% -m %LR% -a %alutfn% -v %hlutfn% -b 2

call %prog% -r %root% -i %input% -m %TT% -T %interval% -I %maxiter% -L %lambda% -S %stop% -b 2
::call %prog% -r %root% -i %input% -m %LL% -a %alutfn% -v %hlutfn% -b 2	



