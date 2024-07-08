#!/usr/bin/env python
## By Shuai Wang
## Notes: Do not tune parmeters, instead use L2 linear-SVM as a reference (Varoquaux et al., 2017)


## Set environment (packages, functions, working path etc.)
# Load up packages
import os
import pickle
import nilearn.decoding
import pandas as pd
import numpy as np
from datetime import datetime
from nilearn.image import load_img, index_img, mean_img, new_img_like
from nilearn.input_data import NiftiMasker
from sklearn import svm
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import LeaveOneGroupOut, PredefinedSplit, cross_validate, permutation_test_score
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_classif
# Setup path
dir_main = '/scratch/swang/agora/CP00/AudioVisAsso'  # the project main folder
dir_data = os.path.join(dir_main, 'derivatives')     # different analyses after pre-processing
dir_mvpa = os.path.join(dir_data, 'multivariate')    # multivariate analyses folder
dir_mask = os.path.join(dir_data, 'masks')           # masks folder
# Read subjects information
f_subjs = os.path.join(dir_mvpa, 'participants_final.tsv')     # subjects info
subjects = pd.read_table(f_subjs).set_index('participant_id')  # the list of subjects
n = len(subjects)
# Task parameters
np.random.seed(21)  # random seed for sklearn
task = 'task-AudioVisAssos1word'    # task name
spac = 'space-MNI152NLin2009cAsym'  # anatomical template that used for preprocessing by fMRIPrep 
mods = ['V', 'A']       # stimulus modalities
## ---------------------------


## MVPA parameters
# Prepare classifiers without tuning parameters
clf_models = [svm.SVC(kernel='linear', max_iter=-1),
              svm.SVC(max_iter=-1),
              LinearDiscriminantAnalysis(),
              GradientBoostingClassifier()]
clf_tokens = ['SVClin', 'SVCrbf', 'LDA', 'GBC']  # classifier abbreviations
nmodels = len(clf_tokens)
# Searchlight parameters
R      = 4   # 57 voxels
N_JOBs = -1  # -1 means all CPUs
# Read ROIs information
f_roi = os.path.join(dir_mvpa, 'group_masks_searchlight.csv')  # ROIs info
df_roi = pd.read_csv(f_roi).set_index('label')              # the list of ROIs
df_roi = df_roi[df_roi.input == 1]                             # only inculde new ROIs
nroi = len(df_roi)
## ---------------------------

now = datetime.now()
print("========== START JOB : %s ==========\n" % now.strftime("%Y-%m-%d %H:%M:%S"))

## MVPA
# Read trial labels
f_labels = os.path.join(dir_mvpa, f'group_{task}_labels-group-trial.tsv')  # trial-wise labels
df_labels = pd.read_table(f_labels) # trial labels
# Do MVPA with leave-one-subject-out CV
f_betas = "%s/group_LSS_nilearn.nii.gz" % dir_mvpa
for imod in mods:
    # Loop ROI
    for iroi in range(nroi):
        thisroi = df_roi.index[iroi]
        f_mask = os.path.join(dir_mask, 'group', "group_%s_mask-%s.nii.gz" % (spac, thisroi))            
        img_mask = load_img(f_mask)
        # Do MVPA with each classifier
        for clf_token, clf_model in zip(clf_tokens, clf_models):
            # Loop participant
            for i in range(n):
                subj = subjects.index[i]
                # Unimodal searchlight
                f_map = os.path.join(dir_mvpa, 'groupMVPA', 'SearchlightMaps', f'{subj}_gMVPA-{clf_token}_LOSOCV_ACC-{imod}_searchlight-{R}mm_mask-{thisroi}.nii.gz')
                if not os.path.exists(f_map):
                    print(f'Searchlight using the classifer {clf_token} within {thisroi} for the modality {imod} for participant {subj}:')
                    labels_uni = np.zeros(len(df_labels))  # initial labels for unimodal decoding
                    labels_train = df_labels['correct'] * (df_labels['modality'] == imod) * (df_labels['participant_id'] != subj)
                    labels_valid = df_labels['correct'] * (df_labels['modality'] == imod) * (df_labels['participant_id'] == subj)
                    labels_all = labels_train | labels_valid
                    labels_uni[labels_train] = -1  # these participants are for unimodal training (-1)
                    CV_uni = PredefinedSplit(labels_uni[labels_all])  # pre-defined CV for unimodal decoding
                    targets = df_labels['lexicon'][labels_all].values   
                    betas = index_img(f_betas, labels_all)           # select betas with this modality
                    betas_mean = mean_img(betas)        # as a template to output results
                    searchlight_uni = nilearn.decoding.SearchLight(mask_img=img_mask, radius=R, estimator=clf_model, n_jobs=N_JOBs, cv=CV_uni)
                    searchlight_uni.fit(betas, targets)
                    searchlight_uni_map = new_img_like(betas_mean, searchlight_uni.scores_)
                    searchlight_uni_map.to_filename(f_map)
                # Cross-modal searchlight
                f_map = os.path.join(dir_mvpa, 'groupMVPA', 'SearchlightMaps', f'{subj}_gMVPA-{clf_token}_LOSOCV_ACC-{imod}2_searchlight-{R}mm_mask-{thisroi}.nii.gz')
                if not os.path.exists(f_map):
                    print(f'Searchlight using the classifer {clf_token} within {thisroi} for the modality {imod}2 for participant {subj}:')
                    labels_cross = np.zeros(len(df_labels))  # initial labels for cross-modal decoding
                    labels_train = df_labels['correct'] * (df_labels['modality'] == imod) * (df_labels['participant_id'] != subj)
                    labels_valid = df_labels['correct'] * (df_labels['modality'] != imod) * (df_labels['participant_id'] == subj)
                    labels_all = labels_train | labels_valid
                    labels_cross[labels_train] = -1  # these are for cross-modal training (-1)
                    CV_cross = PredefinedSplit(labels_cross[labels_all])  # pre-defined CV for cross-modal decoding
                    targets = df_labels['lexicon'][labels_all].values   
                    betas = index_img(f_betas, labels_all)           # select betas with this modality
                    betas_mean = mean_img(betas)        # as a template to output results
                    searchlight_cross = nilearn.decoding.SearchLight(mask_img=img_mask, radius=R, estimator=clf_model, n_jobs=N_JOBs, cv=CV_cross)
                    searchlight_cross.fit(betas, targets)
                    searchlight_cross_map = new_img_like(betas_mean, searchlight_cross.scores_)
                    searchlight_cross_map.to_filename(f_map)

print('Finish Searchlight group MVPA.\n')
## ---------------------------

now=datetime.now()
print("========== ALL DONE : %s ==========\n" % now.strftime("%Y-%m-%d %H:%M:%S"))
