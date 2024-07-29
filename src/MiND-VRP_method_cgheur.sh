#!/bin/bash

#SBATCH -a 1-48
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --constraint=centos7
#SBATCH --mem=32G
#SBATCH -p sched_mit_sloan_batch
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bmair@mit.edu
#SBATCH --output=batch_Met_BHCG\%a.log

OUTPUT_DIR="/home/bmair/DeviatedFixedRouteMicrotransit/data/output"

# microtransit method variant (Benders + Heur CG)
module load julia/1.5.2
module load gurobi/8.1.1
julia MiND-VRP_method_cgheur.jl -o "${OUTPUT_DIR}" -f "trials_MiND-VRP.csv" -t "${SLURM_ARRAY_TASK_ID}"


mv "batch_Met_BHCG${SLURM_ARRAY_TASK_ID}.log" "${OUTPUT_DIR}"
