# Run from the current directory and with current environment
#$ -cwd -V

# Ask for some time (hh:mm:ss max of 48:00:00)
#$ -l h_rt=48:00:00

# Ask for one node
#$ nodes=1

# Send emails when job starts and ends
#$ -m be

# Now run the job
module load R
R CMD BATCH de_analysis.R

