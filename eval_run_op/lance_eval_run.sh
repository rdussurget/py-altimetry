#!/bin/csh
#set echo

source /usr/share/modules/init/csh
module purge
# -- CDAT "home13" PREVIMER
# module use /home13/caparmor/previmer/op/bin/vacumm/modulesfiles
# module load cdat
# module load gcc
# -- CDAT "old"
module use /home11/caparmor/mars/PYTHON/modulefiles
module load cdat
module load gcc
# -- CDAT "git"
# module use /home11/caparmor/mars/PYTHON/modulefiles
# module load cdat-git


pwd
cd /home1/caparmor/sskrypni/pyvalid/trunk/eval_run_op

\rm  eval_run.o*

# update date interval for diag.
./update_config.cfg.py

# launch diag
# 
#en interactif
#./eval_run_sst_op.py
./init_eval_run_sst_op_previmer.py

cd /home/caparmor-work/sskrypni/data_test/SST/SST_SEVIRI
\rm *

cd /home/caparmor-work/sskrypni/data_test/MODEL/MARS_MANGA
\rm *
