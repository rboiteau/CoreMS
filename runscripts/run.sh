#!/bin/bash

for i in 1 2 3 4 5 6 7
do
   RUNFILE=230301_spring-env_pos_21T_assign_0$i.py
   /home/dewey/CoreMS/env/bin/python $RUNFILE
done

python 230301_spring-env_pos_21T_post-assign-proc.py