#!/bin/bash
## By Shuai Wang

# Platform
platform='mesoc'
case "$platform" in
  mesoc)
    mdir='/CP00'                       # the project Main folder @mesocentre
    export PATH="$mdir/nitools:$PATH"  # setup tools if @mesocentre
    njob=16
    ;;
  totti)
    mdir='/data/mesocentre/data/agora/CP00'  # the project Main folder @totti
    njob=4
    ;;
  *)
    echo -e "Please input a valid platform!"
    exit 1
esac
# Setup path
ddir="$mdir/AudioVisAsso"                # experiment Data folder (BIDS put into fMRIPrep)
adir="$ddir/derivatives/afni"            # AFNI output folder
# Processing parameters
readarray subjects < $mdir/CP00_subjects.txt
task='task-AudioVisAssos1word'        # task name
spac='space-MNI152NLin2009cAsym'      # anatomical template that used for preprocessing by fMRIPrep
bold='desc-preproc_bold'              # the token for the preprocessed BOLD data (without smoothing)
regs='desc-confounds_timeseries'      # the token for fMRIPrep output nuisance regressors
anat='desc-preproc_T1w_brain'         # skull-stripped anatomical image
deno='NR24a'                          # denoising strategy
hmpv="dfile_motion_${deno}"           # all head motion NRs that should be regressed out
ortv="dfile_signal_${deno}"           # all non-motion NRs that should be regressed out
cenv='dfile_censor_FD'                # censors
tr=1.148                              # TR in seconds
nrun=5                                # number of runs
ntp=1000                              # total number of TRs
lrun=(200 200 200 200 200)            # number of TRs for each run
drun=(229.6 229.6 229.6 229.6 229.6)  # duration for each run
cons=("WA" "WV" "PA" "PV")            # conditions

# Run LSS for each subject
for subj in ${subjects[@]};do
  echo -e "run LSS for $task for subject : $subj ......"

  wdir="$adir/$subj/$task"                                     # the Working folder
  oglm="${subj}_${task}_GLM.wBIM.wPSC.w${deno}"                # the token for the Output GLM, psc means "percent signal change"
  ldir="$wdir/$oglm/trial-wise_estimates"                      # the LSS working folder
  mask="$wdir/${subj}_${task}_${spac}_desc-brain_mask.nii.gz"  # brain mask that used for 3dLSS
  
  # Create a working folder for LSS
  if [ ! -d $ldir ];then
    mkdir -p $ldir
    cp -r $wdir/$oglm/stimuli $ldir
  fi
  echo -e "# enter the LSS working directory..."
  cd $ldir  # enter the LSS working directory

  # Concatenate all scaled runs
  3dTcat -prefix LSS.${subj}_${task}.all.scale $wdir/$oglm/pb02.${subj}_${task}.r*.scale+tlrc.HEAD  
  
  # Prepare nuisance regressors
  tar -mvxf $wdir/confounds/${oglm}.1D.tar.gz --strip-components 7 -C $wdir/confounds  # use --strip-components 7 to remove folder structure (a bad but only solution)
  cp -r $wdir/confounds/${subj}_${task}_${cenv}.1D $ldir
  1d_tool.py -infile $wdir/confounds/${subj}_${task}_${hmpv}.1D -set_run_lengths ${lrun[*]} -demean -write ${subj}_${task}_${hmpv}_demean.1D
  1d_tool.py -infile ${subj}_${task}_${hmpv}_demean.1D -set_run_lengths ${lrun[*]} -split_into_pad_runs ${subj}_${task}_${hmpv}_demean_all
  for irun in $(seq 1 $nrun);do
    frun=`printf "%02d" $irun`
    1d_tool.py -infile $wdir/confounds/${subj}_${task}_run-${frun}_${ortv}.1D -demean -write ${subj}_${task}_run-${frun}_${ortv}_demean.1D
    1d_tool.py -infile ${subj}_${task}_run-${frun}_${ortv}_demean.1D -pad_into_many_runs $irun $nrun -set_run_lengths ${lrun[*]} -write ${subj}_${task}_run-${frun}_${ortv}_demean_all.1D
  done
  rm -r $wdir/confounds/*.1D
  
  # Generate the global timing
  for icond in stimuli/${subj}_${task}_events-cond?.txt;do
    timing_tool.py -timing $icond -local_to_global ${icond::-4}_global.txt -run_len ${drun[*]}
  done
 
  # Create the design matrix for each condition
  # WA: Words Auditory
  cprefix='con1_WA'
  3dDeconvolve -nodata $ntp $tr \
      -censor ${subj}_${task}_${cenv}.1D \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r01.1D head_motion_run1 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r02.1D head_motion_run2 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r03.1D head_motion_run3 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r04.1D head_motion_run4 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r05.1D head_motion_run5 \
      -ortvec ${subj}_${task}_run-01_${ortv}_demean_all.1D nuisance_regressors_run1 \
      -ortvec ${subj}_${task}_run-02_${ortv}_demean_all.1D nuisance_regressors_run2 \
      -ortvec ${subj}_${task}_run-03_${ortv}_demean_all.1D nuisance_regressors_run3 \
      -ortvec ${subj}_${task}_run-04_${ortv}_demean_all.1D nuisance_regressors_run4 \
      -ortvec ${subj}_${task}_run-05_${ortv}_demean_all.1D nuisance_regressors_run5 \
      -polort 2 \
      -global_times \
      -num_stimts 4 \
      -stim_times_IM 1 stimuli/${subj}_${task}_events-cond1_global.txt GAM -stim_label 1 WA \
      -stim_times 2 stimuli/${subj}_${task}_events-cond2_global.txt GAM -stim_label 2 WV \
      -stim_times 3 stimuli/${subj}_${task}_events-cond3_global.txt GAM -stim_label 3 PA \
      -stim_times 4 stimuli/${subj}_${task}_events-cond4_global.txt GAM -stim_label 4 PV \
      -jobs $njob \
      -x1D ${cprefix}_X.xmat.1D -xjpeg ${cprefix}_X.jpg -x1D_uncensored ${cprefix}_X.nocensor.xmat.1D \
      -x1D_stop
  # WV: Words Visual
  cprefix='con2_WV'
  3dDeconvolve -nodata $ntp $tr \
      -censor ${subj}_${task}_${cenv}.1D \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r01.1D head_motion_run1 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r02.1D head_motion_run2 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r03.1D head_motion_run3 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r04.1D head_motion_run4 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r05.1D head_motion_run5 \
      -ortvec ${subj}_${task}_run-01_${ortv}_demean_all.1D nuisance_regressors_run1 \
      -ortvec ${subj}_${task}_run-02_${ortv}_demean_all.1D nuisance_regressors_run2 \
      -ortvec ${subj}_${task}_run-03_${ortv}_demean_all.1D nuisance_regressors_run3 \
      -ortvec ${subj}_${task}_run-04_${ortv}_demean_all.1D nuisance_regressors_run4 \
      -ortvec ${subj}_${task}_run-05_${ortv}_demean_all.1D nuisance_regressors_run5 \
      -polort 2 \
      -global_times \
      -num_stimts 4 \
      -stim_times 1 stimuli/${subj}_${task}_events-cond1_global.txt GAM -stim_label 1 WA \
      -stim_times_IM 2 stimuli/${subj}_${task}_events-cond2_global.txt GAM -stim_label 2 WV \
      -stim_times 3 stimuli/${subj}_${task}_events-cond3_global.txt GAM -stim_label 3 PA \
      -stim_times 4 stimuli/${subj}_${task}_events-cond4_global.txt GAM -stim_label 4 PV \
      -jobs $njob \
      -x1D ${cprefix}_X.xmat.1D -xjpeg ${cprefix}_X.jpg -x1D_uncensored ${cprefix}_X.nocensor.xmat.1D \
      -x1D_stop
  # PA: Pseudowords Auditory
  cprefix='con3_PA'
  3dDeconvolve -nodata $ntp $tr \
      -censor ${subj}_${task}_${cenv}.1D \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r01.1D head_motion_run1 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r02.1D head_motion_run2 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r03.1D head_motion_run3 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r04.1D head_motion_run4 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r05.1D head_motion_run5 \
      -ortvec ${subj}_${task}_run-01_${ortv}_demean_all.1D nuisance_regressors_run1 \
      -ortvec ${subj}_${task}_run-02_${ortv}_demean_all.1D nuisance_regressors_run2 \
      -ortvec ${subj}_${task}_run-03_${ortv}_demean_all.1D nuisance_regressors_run3 \
      -ortvec ${subj}_${task}_run-04_${ortv}_demean_all.1D nuisance_regressors_run4 \
      -ortvec ${subj}_${task}_run-05_${ortv}_demean_all.1D nuisance_regressors_run5 \
      -polort 2 \
      -global_times \
      -num_stimts 4 \
      -stim_times 1 stimuli/${subj}_${task}_events-cond1_global.txt GAM -stim_label 1 WA \
      -stim_times 2 stimuli/${subj}_${task}_events-cond2_global.txt GAM -stim_label 2 WV \
      -stim_times_IM 3 stimuli/${subj}_${task}_events-cond3_global.txt GAM -stim_label 3 PA \
      -stim_times 4 stimuli/${subj}_${task}_events-cond4_global.txt GAM -stim_label 4 PV \
      -jobs $njob \
      -x1D ${cprefix}_X.xmat.1D -xjpeg ${cprefix}_X.jpg -x1D_uncensored ${cprefix}_X.nocensor.xmat.1D \
      -x1D_stop
  # PV: Pseudowords Visual
  cprefix='con4_PV'
  3dDeconvolve -nodata $ntp $tr \
      -censor ${subj}_${task}_${cenv}.1D \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r01.1D head_motion_run1 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r02.1D head_motion_run2 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r03.1D head_motion_run3 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r04.1D head_motion_run4 \
      -ortvec ${subj}_${task}_${hmpv}_demean_all.r05.1D head_motion_run5 \
      -ortvec ${subj}_${task}_run-01_${ortv}_demean_all.1D nuisance_regressors_run1 \
      -ortvec ${subj}_${task}_run-02_${ortv}_demean_all.1D nuisance_regressors_run2 \
      -ortvec ${subj}_${task}_run-03_${ortv}_demean_all.1D nuisance_regressors_run3 \
      -ortvec ${subj}_${task}_run-04_${ortv}_demean_all.1D nuisance_regressors_run4 \
      -ortvec ${subj}_${task}_run-05_${ortv}_demean_all.1D nuisance_regressors_run5 \
      -polort 2 \
      -global_times \
      -num_stimts 4 \
      -stim_times 1 stimuli/${subj}_${task}_events-cond1_global.txt GAM -stim_label 1 WA \
      -stim_times 2 stimuli/${subj}_${task}_events-cond2_global.txt GAM -stim_label 2 WV \
      -stim_times 3 stimuli/${subj}_${task}_events-cond3_global.txt GAM -stim_label 3 PA \
      -stim_times_IM 4 stimuli/${subj}_${task}_events-cond4_global.txt GAM -stim_label 4 PV \
      -jobs $njob \
      -x1D ${cprefix}_X.xmat.1D -xjpeg ${cprefix}_X.jpg -x1D_uncensored ${cprefix}_X.nocensor.xmat.1D \
      -x1D_stop
  
  # run 3dLSS to obtain trial-wise estimates
  i=1
  for icon in ${cons[@]};do
    echo -e "run LSS for condition $i $icon ......"
    cprefix=`printf "con%d_%s" $i $icon`
    3dLSS -matrix ${cprefix}_X.xmat.1D -input LSS.${subj}_${task}.all.scale+tlrc.HEAD -mask $mask -save1D LSS_${cprefix}.1D -prefix LSS.stats.${subj}_${task}_${cprefix}
    let i+=1
  done

  echo -e "======== everything has been done for subject: $subj ========"
done
