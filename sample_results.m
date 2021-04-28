% 1. Illustration of how regular matrices can be recovered from the sparse
% representation
% 2. Sanity check of the permutation matrices using visualization
% 3. An example of analysis:
%   compare the permutation results on test-retest pairs and on unrelated
%   pairs

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


%% An example of analysis:
% Compare the permutation results on test-retest pairs and on unrelated
% pairs.
% We need to have test and retest results under output/matching_results/
% in order to run the following codes.

% load how retest subjects are matched with test subjects
load('output/test-retest_match_id.mat');

% first retest matrix corresponds to the 21st matrix in the test data
loc(1)

% compare the sum of permutation matrices between (retest 1st - test 21st)
retest_1 = load(strcat('output/matching_results/P_retest_1_', penalty, '.mat'));
sum_swaps = retest_1.sum_swaps;
P_retest_1 = full(sparse(sum_swaps(:, 1), sum_swaps(:, 2), sum_swaps(:, 3), p, p));

test_21 = load(strcat('output/matching_results/P_retest_21_', penalty, '.mat'));
sum_swaps = test_21.sum_swaps;
P_test_21 = full(sparse(sum_swaps(:, 1), sum_swaps(:, 2), sum_swaps(:, 3), p, p));

% difference between retest 1st and test 21st
matched_pair_diff = norm(P_retest_1 - P_test_21, 'fro')

% compare the sum of permutation matrices between (retest 1st - retest 2nd-41st)
mismatched_pair_diff = zeros(1, 997);
for i = 1:997
    test_i = load(strcat('output/matching_results/P_test_', num2str(i) ,'_', penalty, '.mat'));
    sum_swaps = test_i.sum_swaps;
    P_test_i = full(sparse(sum_swaps(:, 1), sum_swaps(:, 2), sum_swaps(:, 3), p, p));
    mismatched_pair_diff(1, i) = norm(P_retest_1 - P_test_i, 'fro');
end

mean(mismatched_pair_diff)
std(mismatched_pair_diff)

% diff(retest_1, test_21) is smaller than the average of diff(retest_1,
% test_i), i=1,...,997
% This is an example of how test-retest pair tends to have permutaiton
% matrices that are more similar than the unrelated pair.