# Run from the current directory and with current environment
#$ -cwd -V

# Ask for some time (hh:mm:ss max of 48:00:00)
#$ -l h_rt=48:00:00

<<<<<<< HEAD
# Ask for one node
#$ nodes=1
=======
# Ask for some memory (by default, 1G, without a request)
#$ -l h_vmem=6G
>>>>>>> 802c9dfb41e4ab844a4bc1b553dcadc5eaf92a15

# Send emails when job starts and ends
#$ -m be

# Now run the job
module load R
R CMD BATCH de_analysis.R

