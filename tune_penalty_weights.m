% Tune penalty weights
% Use test-retest data to tune the penalty weights for yeo_and_rp penalty
% Basically, do a coarse grid search and identify the smallest weight that
% suppresses all swaps
% 
% Chang Su, c.su@yale.edu

% load test-retest data
retest = load('/Users/chang/Documents/research/brain_connectivity/data/precision/test_retest/cc400_regFCprec_concat_hpf_retest41.mat');
test = load('/Users/chang/Documents/research/brain_connectivity/data/precision/cc400_regFCprec_concat_hpf_997subj.mat'); 
[idx,loc] = ismember(retest.subj_id, test.subj);

% set which penalty to look at
penalty = 'yeo_all_pairs'; % penalize for yeo 3 networks and any region-pairs with at least one region invovled in yeo 3 networks

if strcmp(penalty, 'yeo')
    penalty_m = load('output/yeo_penalty.mat').yeo_p;
elseif strcmp(penalty, 'yeo_and_rp')
    penalty_m = load('output/yeo_and_rp_penalty.mat').yeo_and_rp;
elseif strcmp(penalty, 'yeo_and_r')
    penalty_m = load('output/yeo_and_r_penalty.mat').yeo_and_r;
elseif strcmp(penalty, 'yeo_all_pairs')
    penalty_m = load('output/yeo_all_pairs_penalty.mat').yeo_all_pairs;
end    
    
n_iter = 10;

lambdas = (1:1:5)*1e-4; 
nl = 5;
n_off_diag_list = zeros(1, nl);
diff_m_list = zeros(1, nl);
n_penalized_match_list = zeros(1, nl);

for j = 1:nl
    lambda = lambdas(j)
    diff_m = zeros(1, 41);
    P_same = cell(1, 41);
    off_diag_matches = zeros(1, 41);
    penalized_matches = zeros(1, 41);
    
    for i = 1:41
        m1 = squeeze(test.C(loc(i), :, :));
        m1(logical(eye(392))) = 0;
        m2 = squeeze(retest.C(i, :, :));
        m2(logical(eye(392))) = 0;
        [P, diff, obj, new_m, cost_m, sum_swaps] = iterative_procedure(m1, m2, n_iter, lambda, penalty_m);
        P_same{1,i} = P;
        diff_m(i) = diff(max(find(obj)));
        off_diag_matches(1, i) = 392 - sum(P(logical(eye(392))));
        penalized_matches(1, i) = sum(sum(P.*penalty_m));
    end

    % difference between permuted graph and the target graph
    diff_m_list(1, j) = mean(diff_m);
    % number of all off-diagonal swaps
    n_off_diag_list(1, j) = sum(off_diag_matches);
    % number of swaps in the penalized region
    n_penalized_match_list(1, j) = sum(penalized_matches);
end

% check under which penalty did all swaps get muted
diff_m_list
n_off_diag_list
n_penalized_match_list

% visualize the sum of 41 test-retest permutation matrices
% after applying the graph matching algorithm

% note that lambda = 5e-4 gives the same result as lambda = 3e-4,
% which is the optimal weight

P_sum = zeros(392);
for i = 1:41
    P_sum = P_sum + P_same{1,i};
end

yeo_mapping = load('output/yeo_index.mat');

addpath('visualization_tools')
redbluemap = brewermap(100, 'RdBu');
redbluemap = flipud(redbluemap);
bluemap = brewermap(100, 'Blues');
bluemap = flipud(bluemap);
redmap =  brewermap(100, 'Reds');

plot_heatmap(P_sum(yeo_mapping.re_index, yeo_mapping.re_index), [0 1], redmap, ...
    strcat('Sum of test-retest permutation matrices with ', penalty, ' penalty'), yeo_mapping.cluster_count, true(1))

save(strcat('output/matching_results/test_retest_', penalty ,'_sum_swaps.mat'), 'P_sum');
