#!/bin/bash
## By Shuai Wang

# Setup parameters for slurm
#SBATCH -J fMRIPrep
#SBATCH -p skylake
#SBATCH --nodes=1
#SBATCH -A b222
#SBATCH -t 6-12
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH -o ./ps01_PREPROC_fmriprep_%j.out
#SBATCH -e ./ps01_PREPROC_fmriprep_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ws1011001@gmail.com
#SBATCH --mail-type=BEGIN,END
# Setup working path
dir_main='/scratch/swang/agora/CP00'  # the project Main folder
simg='/scratch/swang/simages'     # singularity images folder
# Subjects list
declare -a subjects=$(seq -w 1 22)

# Do preprocessing for each subject
for sid in ${subjects[@]};do
  echo -e “Pre-processing subject: sub-${sid} using fMRIPrep on: $SLURM_NODELIST”
  # Remove temporary files if exist
  if [ "$(ls -A $dir_main/SWAP)" ];then rm -r $dir_main/SWAP/*;fi
  # Preprocessing using fMRIPrep 
  singularity run --cleanenv -B $dir_main:/work $simg/fmriprep-20.2.0 --fs-license-file /work/license.txt \
    /work/AudioVisAsso /work/AudioVisAsso/derivatives participant \
    --participant-label $sid \
    -w /work/SWAP \
    --ignore slicetiming \
    --use-syn-sdc \
    --output-spaces MNI152NLin2009cAsym \
    --fs-no-reconall \
    --return-all-components \
    --stop-on-first-crash \
    --skip_bids_validation
  echo -e "Finish pre-processing for subject :sub-${sid}. Please check it out."
done
# Double-check to clean temporary files
if [ "$(ls -A $dir_main/SWAP)" ];then rm -r $dir_main/SWAP/*;fi

