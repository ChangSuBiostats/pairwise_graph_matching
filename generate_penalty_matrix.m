% Generate penalty matrices
% with the following options:
%   1. vanilla (equal penalty for any swaps)
%   2. hclust two regions only (penalty for four blocks between these two
%   regions)
%   hierarchical cluters)
%
% Chang Su, c.su@yale.edu

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

% Create two regions only penalty matrix
addpath /Users/chang/Documents/research/brain_connectivity/code
hc_func = load('hclust_11_functional.mat');
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