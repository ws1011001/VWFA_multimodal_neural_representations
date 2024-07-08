%% By Shuai Wang

%% Clean up
close all
clear
clc
%% ---------------------------

%% Set environment (packages, functions, working path etc.)
% Initialize brainsotrm
brainstorm nogui   % starts the interface but hides it.
%brainstorm server  % run on headless servers (computation clusters with no screen attached)
% Setup working path
protocol = 'brainstorm_SEEG_LectureVWFA';   % the name of the database
dir_main = '/media/wang/BON/Projects/CP01';     % the project folder
dir_bids = fullfile(dir_main, 'SEEG_LectureVWFA');  % the main working folder 
dir_dev = fullfile(dir_bids, 'derivatives');       % path to the BIDS derivatives
dir_bst = fullfile(dir_dev, protocol);            % path to the BST database
dir_data = fullfile(dir_bst, 'data');              % path to the functional database
% Read the subjects list and information
fid = fopen(fullfile(dir_main, 'subjects_sEEG.txt'));
subjects = textscan(fid, '%s');
fclose(fid);
subjects = subjects{1};  % the subjects list
n = length(subjects);
subjects_info = readtable(fullfile(dir_bids, 'participants.tsv'), 'FileType', 'text');
% Processing parameters
ptoken = 'ses-01_task-RS_run-01_ieeg';     % raw data token
notch_filters = [50, 100, 150, 200, 250];  % notch frequencies in Hz
freq_highpass = 0.3;                       % 0.3 Hz recommended by AnneSo
freq_lowpass = 0;                          % 0 means disable
events.repetition = 'AAp, AAt, AVp, AVt, VAp, VAt, VVp, VVt';   % original eight RSE conditions
events.prime = 'Ap, Vp';                                        % prime conditions merged from the RSE conditions
events_manipulation.merge_prime = [{'AAp, AVp', 'Ap'}; {'VVp, VAp', 'Vp'}];            % merge AAp and AVp as Ap (auditory prime trials), do the same for visual trials
timewindow = [-0.3 0.6];        % epoch window
%% ---------------------------


%% Import raw data
% Create protocol if not exist
if ~exist(dir_bst, 'dir')
  fprintf('Create a new protocol %s in the folder %s. \n', protocol, dir_dev);
  gui_brainstorm('CreateProtocol', protocol, 0, 0, dir_dev);  % create a database in the BIDS derivatives folder
end
% Read data
if ~exist('pFiles', 'var')
  for i = 1:n
    subj = subjects{i};
    sidx = ismember(subjects_info.participant_id, subj);
    f_raw = fullfile(dir_bids, subj, 'ses-01', 'ieeg', sprintf('%s_%s.%s', subj, ptoken, subjects_info.data_format{sidx}));  % raw SEEG recordings (maybe in different formats)
    fprintf('Import raw data file %s into the BST database. \n', f_raw);
    % Read up raw data
    bst_process('CallProcess', 'process_import_data_raw', [], [], 'subjectname', subj, 'datafile', {f_raw, 'SEEG-ALL'}, ...
                'channelreplace', 1, 'channelalign', 0, 'evtmode', 'value');
  end
  % Extract BST filenames for next processes
  pFiles = cellfun(@(x) fullfile(dir_data, x, sprintf('@raw%s_%s/data_0raw_%s_%s.mat', x, ptoken, x, ptoken)), subjects, 'UniformOutput', 0);
end
%% ---------------------------

%% Notch filters : (PSD), notch at 50, 100, 150, 200, 250Hz
if ~exist('sFiles', 'var') || ~isfield(sFiles, 'notch')
  % PSD (to check bad channels)
  bst_process('CallProcess', 'process_psd', pFiles, [], 'timewindow', [], 'win_length', 10, 'win_overlap', 50, ...
              'units', 'physical', 'sensortypes', 'SEEG', 'win_std', 0, ...
              'edit', struct('Comment', 'Power', 'TimeBands', [], 'Freqs', [], 'ClusterFuncTime', 'none', ...
              'Measure', 'power', 'Output', 'all', 'SaveKernel', 0));
  % Notch filters (50, 100, 150, 200, 250Hz)
  sFiles.notch = bst_process('CallProcess', 'process_notch', pFiles, [], 'sensortypes', 'SEEG', 'freqlist', notch_filters, ...
                             'cutoffW', 2, 'useold', 0, 'read_all', 0);
end
%% ---------------------------

%% Band-pass
if ~isfield(sFiles, 'bandpass')
  fprintf('Apply band-pass filtering (%d - %d Hz) on the data. \n', freq_highpass, freq_lowpass);
  % Band-pass filters
  sFiles.bandpass = bst_process('CallProcess', 'process_bandpass', sFiles.notch, [], 'sensortypes', 'SEEG', ...
                                'highpass', freq_highpass, 'lowpass', freq_lowpass, 'tranband', 0, 'attenuation', 'strict', ...  % 60dB
                                'ver', '2019', 'mirror', 0, 'read_all', 0);
end
%% ---------------------------

%% Identify bad channels
% At this stage, examine bad channels and BAD segments based on PSD and continuous data
for i = 1:n
    subj = subjects{i};
    f_dat = bst_process('CallProcess', 'process_select_files_data', [], [], 'subjectname', subj, 'condition', [], 'tag', 'Raw | notch(50Hz 100Hz 150Hz 200Hz 250Hz) | high(0.3Hz)',...
                                 'includebad', 0, 'includeintra', 0, 'includecommon', 0);
    chans_bad = subjects_info.electrodes_bad{i};  % get bad electrodes 
    if ~isempty(chans_bad)  % empty means no bad electrodes
        fprintf('For subject %s, mark the following channels as BAD: %s\n', subj, chans_bad);
        tree_set_channelflag({f_dat.FileName}, 'ClearAllBad');  % cleanup previous records
        tree_set_channelflag({f_dat.FileName}, 'AddBad', chans_bad);
    end    
end
%% ---------------------------

%% Import events
if ~isfield(sFiles, 'events')
  for i = 1:n
    subj = subjects{i};
    fdat = fullfile(dir_data, sFiles.notch(i).FileName);
    fevt = fullfile(dir_bids, subj, 'ses-01', 'ieeg', sprintf('%s_%s_events.csv', subj, ptoken));  % the event file has to be CSV.
    fprintf('Write the final version of events %s into the data %s. \n', fevt, fdat);
    % Write events
    bst_process('CallProcess', 'process_evt_import', fdat, [], 'evtfile', {fevt, 'CSV-TIME'}, 'evtname', '', 'delete', 1);
  end
  sFiles.events = true;  % imported events
end
%% ---------------------------

%% Merge events if needed
events_manip = fieldnames(events_manipulation);  % extract manipulation labels
if ~isfield(sFiles, 'events_created')
  sFiles.events_created = events_manip;
else
  events_manip = events_manip(~ismember(events_manip, sFiles.events_created));  % derive new manips
  sFiles.events_created = [sFiles.events_created; events_manip];                % update sFile
end
if ~isempty(events_manip)
  for i = 1:length(events_manip)
    imanip = events_manip{i};
    evts_manip = events_manipulation.(imanip);
    % Merge trials
    for j = 1:length(evts_manip)
      fprintf('Create a new event %s. \n', evts_manip{j, 2});
      bst_process('CallProcess', 'process_evt_merge', sFiles.bandpass, [], 'evtnames', evts_manip{j, 1}, 'newname', evts_manip{j, 2}, 'delete', 0);
    end
  end
end
%% ---------------------------

%% Epoch
events_epoch = fieldnames(events);
if ~isfield(sFiles, 'trials')
  sFiles.events_epoch = events_epoch;  % update sFiles
else
  events_epoch = events_epoch(~ismember(events_epoch, fieldnames(sFiles.trials)));
  sFiles.events_epoch = [sFiles.events_epoch; events_epoch];
end
if ~isempty(events_epoch)
  % Do epoch for each event
  for i = 1:length(events_epoch)
    ievt = events_epoch{i};
    fprintf('Epoch the event %s. \n', ievt);
    bst_process('CallProcess', 'process_import_data_event', sFiles.bandpass, [], 'subjectname', '', 'condition', '', ...
        'eventname', events.(ievt), 'timewindow', [], 'epochtime', timewindow, 'createcond', 1, 'ignoreshort', 0, 'usectfcomp', 0, 'usessp', 0, 'freq', [], 'baseline', []);
  end
end                           
%% ---------------------------


