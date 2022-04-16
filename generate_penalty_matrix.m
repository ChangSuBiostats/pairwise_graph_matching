% Generate penalty matrices
% with the following options:
%   1. vanilla: equal penalty for any swaps
%   2. hclust two regions only: penalty for four blocks between these two
%   regions
%   3. yeo: penalty for swaps within and between limbic, subcortical and
%   cerrebellum
%   4. yeo_and_rp, yeo and region pairs: in addition to 3, add penalty for any
%   swaps happened in test-retest pairs
%   5. yeo_and_r, yeo and all regions: in addition to 3, add penalty for
%   any regions that have swapped in test-retest pairs
%   6. yeo_all_pairs: in addition to 3, add penalty for any region-pairs
%   that involve regions in 3
% Chang Su, c.su@yale.edu

% load colormaps for visualization
addpath('visualization_tools')
redbluemap = brewermap(100, 'RdBu');
redbluemap = flipud(redbluemap);
bluemap = brewermap(100, 'Blues');
bluemap = flipud(bluemap);
redmap =  brewermap(100, 'Reds');

%% 1. vanilla
% Create vanilla penalty matrix
vanilla_p = zeros(392);
for i = 1:392
    for j = 1:392
        if i ~= j
            vanilla_p(i, j) = 1;
        end
    end
end
save('output/vanilla_penalty.mat', 'vanilla_p');


%% 2. hierarchical clustering based
% Create a matrix to penalize two regions enriched for swappings 
% given by hierarchical clustering
hc_func = load('output/hclust_11_functional.mat');
two_region_p = zeros(392);
for i = 1:392
    for j = 1:392
        if i ~= j
            if any(i == hc_func.left_index(hc_func.T_left == 3)) && ...
                    any(j == hc_func.left_index(hc_func.T_left == 3))
                two_region_p(i, j) = 1;
            end
            if any(i == hc_func.right_index(hc_func.T_right == 1)) && ...
                    any(j == hc_func.right_index(hc_func.T_right == 1))
                two_region_p(i, j) = 1;
            end
            if any(i == hc_func.left_index(hc_func.T_left == 3)) && ...
                    any(j == hc_func.right_index(hc_func.T_right == 1))
                two_region_p(i, j) = 1;
            end
            if any(i == hc_func.right_index(hc_func.T_right == 1)) && ...
                    any(j == hc_func.left_index(hc_func.T_left == 3))
                two_region_p(i, j) = 1;
            end
        end
    end
end

save('output/two_regions_penalty.mat', 'two_region_p');

%% 3. yeo - 3 networks
% Create a matrix to penalize swaps within and between 
% Limbic, Subcortical, Cerebellum under Yeo mapping
yeo_mapping = load('output/yeo_index.mat');

yeo_p = zeros(392);
% limbic: 5, Subcortical: 8, Cerebellum/Brain Stem: 9
three_region_index = ismember(yeo_mapping.network_label, [5, 8, 9]);
yeo_p(three_region_index, three_region_index) = 1;
yeo_p(logical(eye(392))) = 0;

save('output/yeo_penalty.mat', 'yeo_p');

plot_heatmap(yeo_p(yeo_mapping.re_index, yeo_mapping.re_index), [0 1], redmap, ...
    'Yeo penalty', yeo_mapping.cluster_count, true(1));

%% 4. yeo and any region pairs that have swapped in test-retest data
% combine 3 networks and swapped region-pairs
yeo_and_rp = yeo_p;

% load sum of permutation matrices between test-retest pairs
% can be computed with test_retest_swaps.m
P_same_sum = load('output/matching_results/test_retest_no_penalty_sum_swaps.mat').P_same_sum;

yeo_and_rp(P_same_sum > 0) = 1;
yeo_and_rp(logical(eye(392))) = 0;

save('output/yeo_and_rp_penalty.mat', 'yeo_and_rp');

plot_heatmap(yeo_and_rp(yeo_mapping.re_index, yeo_mapping.re_index), [0 1], redmap, ...
    'Yeo_and_rp penalty', yeo_mapping.cluster_count, true(1));


%% 5. yeo and any regions that have swapped in test-retest pairs
n_swaps_per_region = sum(P_same_sum - diag(diag(P_same_sum)), 2);
tab = tabulate(n_swaps_per_region);
bar(tab(:,1), tab(:,3))

yeo_and_r = zeros(392);

% penalize any region pairs swapped more than region_cutoff times to 
region_cutoff = 0;
yeo_and_r(find(n_swaps_per_region >region_cutoff), find(n_swaps_per_region > region_cutoff)) = 1;

% penalize for the three network
yeo_and_r(three_region_index, three_region_index) = 1;
yeo_and_r(logical(eye(392))) = 0;

save('output/yeo_and_r_penalty.mat', 'yeo_and_r');

plot_heatmap(yeo_and_r(yeo_mapping.re_index, yeo_mapping.re_index), [0 1], redmap, ...
    'yeo_and_r', yeo_mapping.cluster_count, true(1));


%% 6. yeo and any region-pairs with at least one region invovled in yeo 3 networks
yeo_all_pairs = zeros(392);
yeo_all_pairs(three_region_index, :) = 1;
yeo_all_pairs(:, three_region_index) = 1;
yeo_all_pairs(logical(eye(392))) = 0;
save('output/yeo_all_pairs_penalty.mat', 'yeo_all_pairs');

plot_heatmap(yeo_all_pairs(yeo_mapping.re_index, yeo_mapping.re_index), [0 1], redmap, ...
    'Yeo 3 networks with all other regions', yeo_mapping.cluster_count, true(1));