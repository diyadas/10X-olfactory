#!/bin/bash
#
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH -t 48:00:00
#SBATCH -p LM
#SBATCH --mem 500GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=diyadas@berkeley.edu,rchance@berkeley.edu

module load gcc

job=$SLURM_JOB_ID
NOW=$(date +"%Y%m%d-%H%M%S")
LOG="test_${job}_${NOW}.Rout"
exec >> "$LOG" 2>&1 || exit 1     # redirect stdout and error to log file, will fail if the logs directory doesn't exist

# as refactored by JW Adams from last commit
# usage: run input_file command to run
run() {
       local input="$1"
       shift;  # remove file from argument list
       echo "command: $* < $input"
       echo "--------------"
       "$@" < "$input"
}

usage() {
       echo "usage: test.sh" >&2
       exit 2
}

[[ $# -eq 0 ]] || usage  # fail if incorrect number of args and print usage info

while true; do free -h >> 'test_'$NOW'_memory.out'; sleep 15; done &

run test.R \
   env R_LIBS=/share/groups/diya-russell/rpack/3.5/ R --vanilla 

