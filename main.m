%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main file to analyse CA-data
%
%   - one class for the stimulus (e.g. chirp)
%   - one class for the CaData
%   - one function file for k-means clustering
%   - one function file for GMM clustering
%   - one function file for functional analysis
%
%   Author: Lars Bollmann, 2019/01/07
%
%   TO Dos: - separate computation & plotting
%           - use "Replicates" option for k-means
%           
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
colordef black;
%--------------------------------------------------------------------------
% Files and dependencies
%--------------------------------------------------------------------------

% T8 superior colliculus 
fname_SC='../01 Data/2018_15_11/T8_181115_0004_stimlog.txt';
SI_Aux_file_SC ='../01 Data/2018_15_11/T8_181115_GC6_00004_SI_Aux_Info.mat';
Ca_file_SC ='../01 Data/2018_15_11/T8_181115_GC6_00004_ExportedCaData.mat';

% T4 visual cortex
fname_VC='../01 Data/2018_15_11/181115_120955_stimlog.txt';
SI_Aux_file_VC ='../01 Data/2018_15_11/T4_181115_GC6_00004_SI_Aux_Info.mat';
Ca_file_VC ='../01 Data/2018_15_11/T4_181115_GC6_ExportedCaData.mat';

addpath('00 functions');
%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Stimulus
%--------------------------------------------------------------------------

% SC
new_chirp_SC = chirp;
new_chirp_SC.readLogFile(fname_SC);
% VC
new_chirp_VC = chirp;
new_chirp_VC.readLogFile(fname_VC);

% return un-scaled default chirp
stim_SC = new_chirp_SC.returnScaledChirp(136);
stim_raw_SC = new_chirp_SC.returnScaledChirp();

%--------------------------------------------------------------------------
% CaData
%--------------------------------------------------------------------------

% SC
data_SC = CaData(SI_Aux_file_SC, Ca_file_SC);
% synchronize stimulus and CaData
data_SC.SyncStimCaData(new_chirp_SC.stim_start_times, new_chirp_SC.stim_end_times)
% filter abnormal cells (arbitrary: mean > 2)
data_SC.cutOffFilter(2)
data_SC.zScoreRows
% calculate average over trials
data_SC.averageOverTrials();
% split trials to get single trial data
data_SC.SplitIntoTrials();
% stack up trials
data_SC.StackUpTrials();
% split trials
data_SC.SplitTrialsIntoParts(3);

% plot data 
%data_SC.plot_all_data_and_stim(stim_raw_SC);

% VC
data_VC = CaData(SI_Aux_file_VC, Ca_file_VC);
% synchronize stimulus and CaData
data_VC.SyncStimCaData(new_chirp_VC.stim_start_times, new_chirp_VC.stim_end_times)
% filter abnormal cells (arbitrary: mean > 2)
data_VC.cutOffFilter(2)
data_VC.zScoreRows
% calculate average over trials
data_VC.averageOverTrials();
% split trials
data_SC.SplitIntoTrials();
% plot data
%data_VC.plot_all_data_and_stim(stim_SC);


%% SNR analysis
%##########################################################################
% calculate SNR for each cell and each trial
data_SC.SNRTrials();
data_SC.plot_SNR_cell_summary(stim_SC);

% calculate SNR for each time point averaged over trials
data_SC.SNRTimePointsCells
data_SC.plot_SNR_per_time_point(stim_raw_SC)


%% k-means clustering
%##########################################################################

% #clusters
nClusters = 10;

% function parameters:
%--------------------------------------------------------------------------
%k_means_clustering(data object,type,nClusters,stimulus,trial selection);

% concatenated trials
%--------------------------------------------------------------------------
%id_k_means_avg = k_means_clustering(data_SC,'temp_avg_trials',nClusters,stim_raw_SC,[]);

% concatenated trials
%--------------------------------------------------------------------------
%id_k_means_temp = k_means_clustering(data_SC,'temp',nClusters,stim_raw_SC,[]);

% cluster localization
%--------------------------------------------------------------------------
%k_means_clustering(data_SC,'temp_cluster_localization',nClusters,stim_raw_SC,[]);

% single trials
%--------------------------------------------------------------------------
%k_means_clustering(data_SC,'temp_single_trial',nClusters,stim_raw_SC,1);

% stacked up trials
%--------------------------------------------------------------------------
%id_k_means_temp_all_trials = k_means_clustering(data_SC,'temp_across_trials',nClusters,stim_raw_SC,[]);

% population states
%--------------------------------------------------------------------------
%k_means_clustering(data_SC,'pop_states',nClusters,stim_raw_SC,[]);

% using split trials
%--------------------------------------------------------------------------
%id_k_means_trial_parts = k_means_clustering(data_SC,"temp_trial_parts",nClusters,stim_raw_SC,[]);

% SNR clustering
%--------------------------------------------------------------------------
%k_means_clustering(data_SC,"SNR_cluster",nClusters,stim_raw_SC,[]);


%% PCA
%##########################################################################
X = data_SC.dFF;
[V, E, D] = pca(X);
D = diag(D);
varEx = cumsum(D) ./ sum(D);

figure;
plot(varEx); 
xlabel('Number of components');
ylabel('Proportion of variance explained');


%% GMM clustering
%##########################################################################

% function parameters:
%--------------------------------------------------------------------------
%GMM_clustering(data matrix,type,nClusters,proportion test set,
% fitting repetitions, array with number of clusters, stimulus)

% #clusters
nClusters = 25;

% proportion of test set
prop_t = 0.2;

% how many times to fit GMM (to minimize influence of initial values)
fitting_rep = 20;

% define what different numbers of clusters should be evaluated
clust_nr_array = 1:3:6;

% input data set
X = data_SC.dFF_all_trials;

% clustering
%--------------------------------------------------------------------------
%id_cl_GMM = GMM_clustering(X,"fitting",nClusters,[],[],[],stim_raw_SC);

% cross-validation
%--------------------------------------------------------------------------
%GMM_clustering(X,"cross_val",[],prop_t,fitting_rep,clust_nr_array,[]);



%% Functional analysis of cells
%##########################################################################

% function parameters:
%--------------------------------------------------------------------------
%FunctionalCellAnalysis(data object,type,id_all_trials,id_avg,id_parts,selected
% average clusters,cell_id_array);

% switching cells for k-means clustering for all trial data
%--------------------------------------------------------------------------
%per_trial_clus = FunctionalCellAnalysis(data_SC,"simple_switching",id_k_means_temp_all_trials,[],[],[]);

% subset of cells by selecting clusters of average data
%--------------------------------------------------------------------------
%FunctionalCellAnalysis(data_SC,"switching_of_subset",id_k_means_temp_all_trials,id_k_means_avg,[],[1,2,14],[]);

% plot switching behaviour of cells
%--------------------------------------------------------------------------
%FunctionalCellAnalysis(data_SC,"plot_switching_cells",id_k_means_temp_all_trials,[],[],[],[1 2]);

% mutual information (cluster,substim) and entropy
%--------------------------------------------------------------------------
%FunctionalCellAnalysis(data_SC,"mutual_info_ent",[],[],id_k_means_trial_parts,[],[]);



%% Functional analysis of clusters
%##########################################################################

% function parameters:
%--------------------------------------------------------------------------
%FunctionalCellAnalysis(data object,type,idx,per trial cluster
%distribution)

% using entropy for trials
%--------------------------------------------------------------------------
%FunctionalClusterAnalysis(data_SC,"char_clusters",id_k_means_temp_all_trials,per_trial_clus);

% transitions between clusters
%--------------------------------------------------------------------------
%FunctionalClusterAnalysis(data_SC,"transitions_clusters",id_k_means_temp_all_trials,per_trial_clus);


