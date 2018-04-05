#!/bin/bash

#$ -S /bin/bash
#$ -q course.q

cd $HOME/Assign_3

#$ -cwd

spades.py --careful -o spades_qsub_T9 -1 /home/nih13008/Assign_3/T9_fwd_paired.fastq -2 /home/nih13008/Assign_3/T9_rvs_paired.fastq
