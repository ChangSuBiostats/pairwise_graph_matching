% 1. Illustration of how regular matrices can be recovered from the sparse
% representation
% 2. Sanity check of the permutation matrices using visualization
%

penalty_m = load('output/yeo_penalty.mat').yeo_p;

i_ref=1;
ref_popu='test';
penalty = 'yeo';
load(strcat('output/matching_results/P_', ref_popu,'_', num2str(i_ref), '_', penalty, '.mat'));

%% recover the permutation matrix from its sparse representation
p = 392;

% recover the permutation matrix for matching subject 1 to subject 10
j = 10;
P_to_10 = full(sparse(swap_positions(j, :, 1),swap_positions(j, :, 2), 1, p, p));

% recover the sum of permutation matrices across matching to subject 1-997
P_sum=full(sparse(sum_swaps(:, 1), sum_swaps(:, 2), sum_swaps(:, 3), p, p));

%% sanity check on the computed results: visualize the permutation matrices

% load ordering of yeo mapping
addpath('visualization_tools')
yeo_mapping = load('output/yeo_index.mat');

% initiate red blue map for plotting FC maps and red maps for plotting
% permutation matrices

redmap =  brewermap(100, 'Reds');

% illustrate plot_heatmap(): plot penalty matrix
penalty_m = load('output/yeo_penalty.mat').yeo_p;
plot_heatmap(penalty_m(yeo_mapping.re_index, yeo_mapping.re_index), [0 1], redmap, ...
    'Yeo penalty', yeo_mapping.cluster_count, true(1));

% permutation matrix of matchin gsubject 1 to subject 10
plot_heatmap(P_to_10(yeo_mapping.re_index, yeo_mapping.re_index), [0 1], redmap, ...
    'subject 1-subject 10', yeo_mapping.cluster_count, true(1));

% sum of permutation matrices
% set the range to be [0,50] swaps
plot_heatmap(P_sum(yeo_mapping.re_index, yeo_mapping.re_index), [0 50], redmap, ...
    'sum of subject 1-all other subjects (max=50)', yeo_mapping.cluster_count, true(1));
% set the range to be [0,1] swaps
% i.e. red color denotes 'at least one swaps when matching to 997 subjects'
plot_heatmap(P_sum(yeo_mapping.re_index, yeo_mapping.re_index), [0 1], redmap, ...
    'sum of subject 1-all other subjects (max=1)', yeo_mapping.cluster_count, true(1));
% It should be clear from this figure that all permutations in the
% limbic-subcortical-cerebellum region are 'muted'


