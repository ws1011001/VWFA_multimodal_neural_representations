#!/bin/bash
## By Shuai Wang

# Platform
platform='mesoc'
case "$platform" in
	mesoc)
		dir_main='/CP00'                       # the project Main folder @mesocentre
  	  	export PATH="$dir_main/nitools:$PATH"  # setup tools if @mesocentre
  	  	njob=16
  	  	;;
  	local)
		dir_main='/data/mesocentre/data/agora/CP00'  # the project Main folder @totti
  	  	njob=4
  	  	;;
  	*)
		echo -e "Please input a valid platform!"
  	  	exit 1
esac
# Setup path
dir_data="$dir_main/AudioVisAsso"                # experiment Data folder (BIDS put into fMRIPrep)
dir_afni="$dir_data/derivatives/afni"            # AFNI output folder
dir_mask="$dir_data/derivatives/masks"
# Processing parameters
readarray subjects < $dir_main/subjects.txt
readarray rois < $dir_afni/ROIs.txt
task='task-AudioVisAssos2words'   # task name
spac='space-MNI152NLin2009cAsym'  # anatomical template that used for preprocessing by fMRIPrep
deno='NR24a'                      # denoising strategy
cons=("SISMa" "SISMv" "SIDMa" "SIDMv" "DISMa" "DISMv" "DIDMa" "DIDMv" "catch")

# Extract PSC
f_psc="$dir_afni/group_${task}_RSE_PSC.csv"
echo "participant_id,ROI_label,condition,PSC" >> $f_psc
for subj in ${subjects[@]};do
	echo -e "extract PSC for $task for subject : $subj ......"
  	dir_task="$dir_afni/$subj/$task"                       # the Working folder
  	oglm="${subj}_${task}_GLM.wBIM.wPSC.w${deno}"  # the token for the Output GLM
  	# Extract PSC for each ROI 
  	f_glm="$dir_task/$oglm/stats.${subj}_${task}+tlrc."
  	for iroi in ${rois[@]};do
		if [ "${iroi::1}" = 'i' ];then
			f_roi="$dir_mask/$subj/${subj}_${spac}_mask-${iroi}.nii.gz"
  	  	else
			f_roi="$dir_mask/group/group_${spac}_mask-${iroi}.nii.gz"
  	  	fi
  	  	i=1
  	  	for icon in ${cons[@]};do
  	  	  	x=$(3dmaskave -q -mask $f_roi ${f_glm}[$i])
  	  	  	echo -e "$subj,$iroi,$icon,$x" >> $f_psc
  	  	  	let i+=3
  	  	done
  	done
done

