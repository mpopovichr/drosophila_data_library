;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.
 
#!/bin/bash -i 
#$ -S /bin/bash 
# 
# MPI-PKS script for job submission script with ’qsub’. 
# Syntax is Bash with special qsub-instructions that begin with ’#$’. 
# For more detailed documentation, see 
#       http://www/closed/getting_started/queuing_system.html 
# 
# (last change of this file: $Id: f20009313e90b5cdb7168260caf8bf57d96b4f26 $) 
 
# --- Mandatory qsub arguments 
# Hardware requirements. 
#$ -l h_rss=2000M,h_fsize=10000M,h_cpu=20:00:00,hw=x86_64 
 
# --- Optional qsub arguments 
# Change working directory - your job will be run from the directory 
# that you call qsub in.  So stdout and stderr will end up there. 
#$ -cwd 
 
# --- Job Execution 
# For faster disk access copy files to /scratch first. 
 
scratch=/scratch/$USER/$$ 
mkdir -p $scratch 
cd $scratch 
cp -r $HOME/Documents/Work/Projects/drosophila_wing_analysis/drosophila_data_library/globals.py . 
 
# Execution - running the actual program. 
# [Remember: Don’t read or write to /home from here.] 
echo "Running␣on␣$(hostname)" 
echo "We␣are␣in␣$(pwd)" 

chmod u+x globals.py
./globals.py
 
# Finish - Copy files back to your home directory, clean up. 
cp -r $scratch $HOME/Documents/Work/Projects/drosophila_wing_analysis/drosophila_data_library/job_data/      # Better use a subdirectory of $HOME. 
cd 
rm -rf $scratch
