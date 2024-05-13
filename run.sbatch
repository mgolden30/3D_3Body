#!/bin/bash
#SBATCH -Jthree_body_search           #Job Name
#SBATCH --account=gts-dpsaltis3
#SBATCH -N1 --ntasks-per-node=1       #number of nodes and cores per node required
#SBATCH --mem-per-cpu=1G              #Memory per core
#SBATCH -t15                          #duration of the job (Ex: 15 mins)
#SBATCH -qembers                      #QOS Name
#SBATCH -oReport-%j.out               #output and error messages
#SBATCH --mail-type=BEGIN,END,FAIL    #Mail preferences
#SBATCH --mail-user=mgolden30@gatech.edu


# copy main executable over from long-term storage
MY_SCRATCH="/storage/home/hcoda1/5/mgolden30/scratch/"
cp /storage/home/hcoda1/5/mgolden30/p-dpsaltis3-0/Matts_data_wonderland/3D_3Body/main  /storage/home/hcoda1/5/mgolden30/scratch/main

cd /storage/home/hcoda1/5/mgolden30/scratch/

chmod u+x main

#module load anaconda3/2-22.05 #Load any modules

srun main
