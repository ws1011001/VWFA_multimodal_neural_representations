%% By Shuai Wang
%% Notes: This script is modified from ASD's script seeg_pipeline_MAIN.m


%% Clean up
close all
clear
clc
%% ---------------------------

%% Set environment (packages, functions, working path etc.)
% % if use HPC
% dir_main = '/CP01';
% addpath(genpath('/CP01/mia'));
% Setup working path
protocol = 'mia_SEEG_LectureVWFA';
dir_main = '/media/wang/BON/Projects/CP01';     % the project folder
dir_bids = fullfile(dir_main, 'SEEG_LectureVWFA');  % the main working folder 
dir_dev  = fullfile(dir_bids, 'derivatives');       % path to the BIDS derivatives
dir_mia  = fullfile(dir_dev, protocol);            % path to the MIA database
% Read the subjects list
fid = fopen(fullfile(dir_main,'subjects_SEEG_lvOT.txt'));
subjects = textscan(fid,'%s');
fclose(fid);
subjects = subjects{1};  % the subjects list
n = length(subjects);
% Processing parameters
freq_step = 10;   % frequency decomposition step (Hz)
freq_bands = [70, 150];
stats.nperm = 1000;  % permutation number for between-condition test
stats.smooth = 20;   % this number x 3 < number of ROIs (for successful smoothing); currently removed
stats.threshp = 0.05;  % p-value for between-condition test
stats.baseline = [-0.3,-0.01];  % baseline from -200ms to -10ms
timewindow = [-0.3,0.6];  % time-window of the epoch
conditions = {'AAp', 'AAt', 'AVp', 'AVt', 'VAp', 'VAt', 'VVp', 'VVt', ...  % RSE conditions
              'Ap', 'Vp'};                                                 % Prime trials to check the contacts bimodal activity
ncond = length(conditions);
pairs = {{'AAp', 'AAt'}, {'AVp', 'AVt'}, {'VAp', 'VAt'}, {'VVp', 'VVt'}};  % RSE contrasts
%% ---------------------------


%% Import data from brainstorm database
% Create protocol if not exist
if ~exist(dir_mia, 'dir')
  fprintf('Create a MIA protocol %s in the folder %s. \n', protocol, dir_dev);
  mkdir(dir_dev, protocol);  % create a database in the BIDS derivatives folder
end
% Read data
bst_db = fullfile(dir_dev, 'brainstorm_SEEG_LectureVWFA', 'data');
for i = 1:n
    subj = subjects{i};
    dir_subj = fullfile(dir_mia, subj);
    if ~exist(dir_subj, 'dir'); mkdir(dir_subj); end
    fprintf('Import BST data into MIA to extract LFP signals for subject %s. \n', subj);
    % Copy BST data to the MIA database
    for j = 1:ncond
      icond = conditions{j};
      if ~exist(fullfile(dir_subj, icond), 'dir')
          mkdir(dir_subj, icond); % create a folder for this condition
          copyfile(fullfile(bst_db, subj, icond), fullfile(dir_subj, icond));
      end
    end
    % Import BST data into MIA
    config.maindir = dir_subj;
    config.outdir = dir_subj;
    mia_s1_extract_bst_data(config);  % edited mia_s1_extract_bst_data
end
%% ---------------------------


%% Extract frequency
for i = 1:n
    subj = subjects{i};
    dir_subj = fullfile(dir_mia, subj);  
    freql = freq_bands(1);  % lower band
    frequ = freq_bands(2);  % upper band           
    fprintf('Extract oscillation amplitudes within frequency band %d to %d Hz for subject %s. \n', freql, frequ, subj);
    % define config
    config.outdir       = dir_subj;
    config.removeEvoked = true;  % remove evoked responses
    config.freqs        = freql:freq_step:frequ;
    config.modetf       = 'Morlet';%'Hilbert'; 
    config.ncycles      = 7;  % Not useful for Hilbert       
    config.mtg = 'bipolar';  % TFA on bipolar signals
    mia_s4_compute_tf_bandwise(config);      
end
%% ---------------------------


%% Single condition statistics to select contacts
for i = 1:n
    subj = subjects{i};
    dir_subj = fullfile(dir_mia, subj);  
    freql = freq_bands(1);  % lower band
    frequ = freq_bands(2);  % upper band           
    for icond = 1:ncond
        dir_cond = fullfile(dir_subj, conditions{icond});
        fprintf('Estimate one-sample significance for frequency band %d to %d Hz for condition %s for subject %s.\n', freql, frequ, conditions{icond}, subj);
        f_tfa = fullfile(dir_cond, sprintf('%s_bipolar_morlet_data_%d_%d_%d_removeEvoked.mat', conditions{icond}, freq_step, freql, frequ));
        mia_s5_compute_stats(f_tfa, stats);  % edited mia_s5_compute_stats     
    end
end
%% ---------------------------


%% Generate ROI table for MIA
% Initialize ROI table
mni_varnames = {'Patient', 'bipolar_label', 'mni_x', 'mni_y', 'mni_z', 'AAL', 'Contact', 'Lateralization', 'Region'};
mni_table    = array2table(zeros(0, length(mni_varnames)), 'VariableNames', mni_varnames);
for i = 1:n
    subj = subjects{i};
    fprintf('Extract ROI information for the subject %s.\n', subj);
    % Read the table of MNI coordinates of bipolar channels
    f_elec = fullfile(dir_bids, subj, 'anat/elecbipolar2atlas.mat');
    d_elec = load(f_elec);
    % Group individual coordinates of contacts
    n_contacts = length(d_elec.coi.label);
    [mia_labels, mia_laterality] = roi_labels4mia(d_elec.coi.label, "A'_1_2");
    mni_cells = [repmat({subj}, n_contacts, 1), ...     % Patient
                 d_elec.coi.label, ...                  % bipolar_label
                 num2cell(d_elec.coi.elecpos_mni), ...  % MNI x, y ,z
                 d_elec.AAL.label(:, 1), ...            % AAL labels
                 mia_labels, ...                        % Contact (MIA contact labels)
                 mia_laterality, ...                    % Lateralization
                 d_elec.AAL.label(:, 1)];               % Region (AAL labels for MIA)
    mni_table = [mni_table; cell2table(mni_cells, 'VariableNames', mni_varnames)];    
end
% Save group data
fmni = fullfile(dir_bids, 'group_space-MNI_contacts-bipolar_ROI-AAL.csv');
writetable(mni_table, fmni);
%% ---------------------------


%% Group between-condition statistics (Multi-patient permutation tests)
% Read ROI table
f_roi = fullfile(dir_mia, 'ROIs_VWFA.mat');
load(f_roi);
% Define parameters
OPTIONS.mia_db = dir_mia;
OPTIONS.dtoken = 'bipolar_morlet_data_10_70_150_removeEvoked';
OPTIONS.win_noedges = [-300, 600];
OPTIONS.timewindow = [-0.3, 0.6];
OPTIONS.nperm = 1000;
OPTIONS.smoth = 20;
OPTIONS.threshp = 0.05;
OPTIONS.Tlim = 10; 
pairs = {{'VVp', 'VVt'}, {'AAp', 'AAt'}};  % RSEs
for ipair = 1:length(pairs)
    cpair = pairs{ipair};
    OPTIONS.contrast = cpair;
    % Extract signals for this ROI
    roi = grp_signals(roi_table, 'Left vOT', OPTIONS);
    % Get baseline from permutation tests
    roi_con = roi_stats_permutations_contrast(roi, OPTIONS);
    % Plot results
    OPTIONS.colors = jet(numel(roi_con.idPt));
    hdur = roi_plot_conditions_fdr2(roi_con, OPTIONS);
end
% figure configuration
cfgBGA.fontsize = 24;
cfgBGA.linewidth = 4;
cfgBGA.smoothing = 1;
cfgBGA.xscale = 0.1000;
cfgBGA.legend = 'none';
cfgBGA.sigbar = true;
cfgBGA.channel = 'HFA';
cfgBGA.axiswidth = 3;
cfgBGA.Time = roi.time;
cfgBGA.TimeLim = 1:length(hdur);  % [-300, 600] ms
% VVp vs. VVt
cfgBGA.YLim = [-0.1 0.6];
cfgBGA.yscale = 0.3000;
signals = [roi_con.sig1_mean, roi_con.sig2_mean]';
ribbons = [roi_con.sig1_ste, roi_con.sig2_ste]';
sig_segment = {[0.183, 0.405]};
plot_timecourse(signals, ribbons, sig_segment, [], 'BGA', cfgBGA, fullfile(dir_mia, 'ROIs_VWFA_VVp-VVt.png'))
% set4a: AAp vs. AAt
cfgBGA.YLim = [-0.1 0.25];
cfgBGA.yscale = 0.1500;
signals = [roi_con.sig1_mean, roi_con.sig2_mean]';
ribbons = [roi_con.sig1_ste, roi_con.sig2_ste]';
sig_segment = {[0.219, 0.254]; [0.392, 0.431]};
plot_timecourse(signals, ribbons, sig_segment, [], 'BGA', cfgBGA, fullfile(dir_mia, 'ROIs_VWFA_AAp-AAt.png'))
%% ---------------------------

