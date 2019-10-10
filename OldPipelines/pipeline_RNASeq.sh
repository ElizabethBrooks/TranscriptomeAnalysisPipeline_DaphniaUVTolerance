#!/bin/bash

qsub -N "job1" foo.sh
qsub -hold_jid "job1" -N "job2" foo.sh