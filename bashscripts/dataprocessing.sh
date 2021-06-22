# Run from the current directory and with current environment
#$ -cwd -V

# Ask for some time (hh:mm:ss max of 48:00:00)
#$ -l h_rt=12:00:00

# Ask for some memory (by default, 1G, without a request)
#$ -l h_vmem=12G

# Ask for output and error files to be sent to a specific folder
#$ -o /nobackup/bs20chlb/outdir/dataprocessing/job-$JOB_ID.stdout
#$ -e /nobackup/bs20chlb/outdir/dataprocessing/job-$JOB_ID.stderr

# Send emails when job starts and ends
#$ -m be

# Now run the job
module load R
R CMD BATCH  /nobackup/bs20chlb/scripts/isoformproject/dataprocessing/Transcript_Filter_Logic/filtering_counts.R /nobackup/bs20chlb/outdir/dataprocessing/R-${JOB_ID}.Rout
