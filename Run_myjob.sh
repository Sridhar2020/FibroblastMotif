
#$ -m be
#$ -M s.seshan@sheffield.ac.uk
#$ -l h_rt=1:00:00
#$ -l mem=4G
#$ -pe openmp 1
#$ -v OMP_NUM_THREADS=1
./aout >> Output

