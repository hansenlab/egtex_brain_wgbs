# Run LDScore estimation.
# Peter Hickey
# 2020-04-20

# SGE variables  ---------------------------------------------------------------

#$ -l h_fsize=500G
#$ -l cegs
#$ -pe local 10

# Load modules -----------------------------------------------------------------

module load conda_R/3.6.x
module load python/2.7.9
export OMP_NUM_THREADS=1

# Run R script -----------------------------------------------------------------

Rscript run_LDScore_estimation.R ${SGE_TASK_ID}
