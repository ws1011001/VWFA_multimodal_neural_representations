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
dir_data="$dir_main/AudioVisAsso"       # experiment Data folder (BIDS put into fMRIPrep)
dir_afni="$dir_data/derivatives/afni"   # AFNI output folder
dir_mask="$dir_data/derivatives/masks"  # masks folder
# Processing parameters
readarray subjects < $dir_main/subjects.txt
readarray rois < $dir_afni/ROIs.txt
task='task-LocaAudio2p5'          # task name
spac='space-MNI152NLin2009cAsym'  # anatomical template that used for preprocessing by fMRIPrep
#mask="$dir_mask/group/group_${spac}_mask-gm0.2_res-${task}.nii.gz"  # GM mask
mask="$dir_mask/group/group_${spac}_mask-gm-lVP_res-${task}.nii.gz"  # left ventral pathway mask
models=("GLM.wPSC.wNR24a")
conditions=("words" "pseudowords" "scrambled")
contrasts=("words-pseudowords" "words-scrambled" "pseudowords-scrambled")  # contrast labels

echo -e "========== START JOB at $(date) =========="

# Extract beta coefficients
for subj in ${subjects[@]};do
	echo -e "extract beta maps for $task for subject : $subj ......"
  	dir_task="$dir_afni/$subj/$task"
  	# Specify stat files
  	for model in ${models[@]};do
		oglm="${subj}_${task}_${model}"
  	  	f_stat="$dir_task/$oglm/stats.${subj}_${task}+tlrc.HEAD"
  	  	# Extract coef maps for group analysis
  	  	i=0
  	  	for cond in ${conditions[@]};do
			this_cond="${cond}#0_Coef"
			f_coef="$dir_task/$oglm/stats.beta_${oglm}_${cond}.nii.gz"
  	  	  	if [ ! -f $f_coef ];then
				3dbucket -fbuc -prefix $f_coef $f_stat[$this_cond]
  	  	  	fi
  	  	  	let i+=1
  	  	done
  	done
done

# Group T-tests
dir_task="$dir_afni/group/$task"
if [ ! -d $dir_task ];then mkdir -p $dir_task;fi
for model in ${models[@]};do
	for cond in ${conditions[@]};do
		# Stack up subjects for group analysis
  	  	f_coef="$dir_task/stats.beta_group_${task}_${model}_${cond}.nii.gz"
  	  	if [ ! -f $f_coef ];then
			3dbucket -fbuc -aglueto $f_coef $dir_afni/sub-*/$task/sub-*_${task}_${model}/stats.beta_sub-*_${cond}.nii.gz
  	  	fi
  	done  
  	# Between-condition T-test with FWE estimation
	for cont in ${contrasts[@]};do
		cond_pair=(${cont//-/ })
		f_con1="$dir_task/stats.beta_group_${task}_${model}_${cond_pair[0]}.nii.gz"
		f_con2="$dir_task/stats.beta_group_${task}_${model}_${cond_pair[1]}.nii.gz"
		f_test="$dir_task/stats.lVP.group_${task}_${model}_TTest_${cont}"
		f_resid="$dir_task/stats.lVP.group.resid_${task}_${model}_TTest_${cont}+tlrc"
		f_acf="$dir_task/stats.lVP.group.ACF_${task}_${model}_TTest_${cont}"
		f_sim="$dir_task/stats.lVP.group.ACFc_${task}_${model}_TTest_${cont}"
		f_fwe="$dir_task/stats.lVP.group.FWE_${task}_${model}_TTest_${cont}"
		if [ ! -f "${f_acf}.1D" ];then
			echo -e "Perform paired T-test for the contrast ${cond_pair[0]} vs. ${cond_pair[1]}."
			# Perform paired T-test
			3dttest++ -setA $f_con1 -setB $f_con2 -mask $mask -exblur 6 -paired -prefix $f_test -resid $f_resid
			# Estimate ACF
			3dFWHMx -ACF -mask $mask -input $f_resid >> ${f_acf}.1D
			mv $dir_main/scripts/3dFWHMx.1D ${f_sim}.1D
			mv $dir_main/scripts/3dFWHMx.1D.png ${f_sim}.png
			# Simulate FWE using 3dClustSim
			read -ra acf <<< $(sed '4!d' ${f_acf}.1D)
			3dClustSim -mask $mask -acf ${acf[0]} ${acf[1]} ${acf[2]} -athr 0.05 0.01 0.005 0.001 -prefix $f_fwe
		fi
	done
done

# Extract coef (PSC) for ROIs
model='GLM.wPSC.wNR24a'
f_psc="$dir_afni/stats.beta_group_${task}_${model}_PSC.csv"
if [ ! -f $f_psc ];then
	echo "participant_id,ROI_label,condition,PSC" >> $f_psc
	for subj in ${subjects[@]};do
		echo -e "Extract beta values (PSC) for subject $subj."
	  	dir_task="$dir_afni/$subj/$task"
	  	oglm="${subj}_${task}_${model}"
	  	# Extract PSC data for each ROI
	  	for iroi in ${rois[@]};do
			if [ "${iroi::1}" = 'i' ];then
				f_roi="$dir_mask/$subj/${subj}_${spac}_mask-${iroi}.nii.gz"
	  	  	else
				f_roi="$dir_mask/group/group_${spac}_mask-${iroi}.nii.gz"
	  	  	fi
			# Resample the ROI to the task resolution
		  	f_roi_tmp="$dir_mask/group/group_${spac}_mask-${iroi}_tmp.nii.gz"
	  	  	3dresample -master $mask -input $f_roi -prefix $f_roi_tmp	
			# Extract PSC for each condition
	  	  	for cond in ${conditions[@]};do
				f_coef="$dir_task/$oglm/stats.beta_${oglm}_${cond}.nii.gz"
				psc=$(3dBrickStat -mean -mask $f_roi_tmp $f_coef)
				echo -e "$subj,$iroi,$cond,$psc" >> $f_psc
	  	  	done
			rm -r $f_roi_tmp
	  	done
	done
fi

echo -e "========== ALL DONE! at $(date) =========="
