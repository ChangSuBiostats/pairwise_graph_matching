%% load data 
retest = load('path_to_retest_data'); % /Users/chang/Documents/research/brain_connectivity/data/precision/test_retest/cc400_regFCprec_concat_hpf_retest41.mat
test = load('path_to_test_data'); % /Users/chang/Documents/research/brain_connectivity/data/precision/cc400_regFCprec_concat_hpf_997subj.mat
[idx,loc] = ismember(retest.subj_id, test.subj);

% load codes for visualization
addpath('visualization_tools')
redbluemap = brewermap(100, 'RdBu');
redbluemap = flipud(redbluemap);
bluemap = brewermap(100, 'Blues');
bluemap = flipud(bluemap);
redmap =  brewermap(100, 'Reds');

%% Run graph matching
n_iter = 20;
lambda = 0; % no regularization
penalty_m = 0;

diff_m = zeros(1, 41);
P_same = cell(1, 41);
obj = zeros(41 ,n_iter+1);

for i = 1:41
    m1 = squeeze(test.C(loc(i), :, :));
    m1(logical(eye(392))) = 0;
    m2 = squeeze(retest.C(i, :, :));
    m2(logical(eye(392))) = 0;
    [P, diff, ob, new_m, cost_m, sum_swaps] = iterative_procedure(m1, m2, n_iter, lambda, penalty_m);
    P_same{1,i} = P;
    diff_m(i) = diff(n_iter + 1);
    obj(i, :) = ob;
end

% summarize 41 test-retest permutation matrices into one swap count matrix
P_same_sum = zeros(392, 392);
for i = 1:41
    P_same_sum = P_same{1, i} + P_same_sum;
end

%% visualize swap count matrix
% load Yeo mapping for re-ordering genes
yeo_mapping = load('output/yeo_index.mat');

% visualize with plot_heatmap under visualization_tools;
plot_heatmap(P_same_sum(yeo_mapping.re_index, yeo_mapping.re_index), [0 4], redmap, ...
    'Test-retest matching (in 4)', yeo_mapping.cluster_count, true(1));

plot_heatmap(P_same_sum(yeo_mapping.re_index, yeo_mapping.re_index), [0 1], redmap, ...
    'Test-retest matching (at least once)', yeo_mapping.cluster_count, true(1));

plot_heatmap(P_same_sum(yeo_mapping.re_index, yeo_mapping.re_index), [0 41], redmap, ...
    'Test-retest matching (in 41)', yeo_mapping.cluster_count, true(1));


%% other exploration
% check number of iterates it takes to converge
conv_iter = zeros(1, 41);
for i = 1:41
    % relative difference
    tmp = find(obj(i,:) == 0);
    conv_iter(i) = tmp(1) - 1;
end
tabulate(conv_iter)