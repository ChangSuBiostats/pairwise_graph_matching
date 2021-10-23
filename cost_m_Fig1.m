% Three example cost matrix with match_to_a_subj 
% used in Figure 1, panel B
% modified based on match_to_a_subj.m
%
% Author:      Chang Su
% email:        c.su@yale.edu
% Matlab ver.:  R2019b
% Date:         23-Oct-2021

i_ref = 1;
lambda = 3e-4;
penalty = 'yeo';
% generated on Chang's laptop
data_path = '/Users/chang/Documents/research/brain_connectivity/data/precision';
ref_popu = 'test';

%% set up paths and directories
% add the path where data are stored
addpath(data_path)
mkdir('output/Fig1')

%% load data
test = load('/Users/chang/Documents/research/brain_connectivity/data/precision/cc400_regFCprec_concat_hpf_997subj.mat');
n = size(test.subj, 2);

%% start matching
    
    % set penalty matrices
    if strcmp(penalty, 'vanilla')
        penalty_m = load('output/vanilla_penalty.mat').vanilla_p;
    elseif strcmp(penalty, 'two_region')
        penalty_m = load('output/two_regions_penalty.mat').two_region_p;
    elseif strcmp(penalty, 'yeo')
        penalty_m = load('output/yeo_penalty.mat').yeo_p;
    end    
    
    % set the FC map from subj i in ref_popu as the reference map
    if strcmp(ref_popu, 'test')
        m_ref = squeeze(test.C(i_ref, :, :));
    elseif strcmp(ref_popu, 'retest')
        m_ref = squeeze(retest.C(i_ref, :, :));
    end
    m_ref(logical(eye(392))) = 0; % remove the effect of diagonal values
    
    % set number of iterations to 1 to save computational time
    n_iter = 1;
    
    % keep track of
    % 1) number of off diagonal swaps
    offdiag_swap_counts = zeros(3,1); ;
    % 2) array of permutation matrices 
    cost_m_array = zeros(392, 392, 3);
    
    inds = [2,3,n];
    for i = 1:3
        m2 = squeeze(test.C(inds(i), :, :));
        m2(logical(eye(392))) = 0;
        [P, d, ~, ~, cost_m] = iterative_procedure(m_ref, m2, n_iter, lambda, penalty_m);
        offdiag_swap_counts(i) = 392 - sum(P(logical(eye(392))));
        cost_m_array(:, :, i) = cost_m;
    end
    
    save('output/Fig1/cost_m_1-2_3_997.mat', 'cost_m_array');