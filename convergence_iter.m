% Plot convergence curve
% for matching subject 1 against all other subjects
% modified from match_to_a_subj.m
%
% Author:      Chang Su
% email:        c.su@yale.edu
% Matlab ver.:  R2019b
% Date:         9-March-2022

addpath '/Users/chang/Documents/research/brain_connectivity/pairwise_matching'

% run the iterative procedure for 10 iterations
n_iter = 20;
i_subj = 1;
match_to_a_subj(i_subj,3e-4,'two_region', '/Users/chang/Documents/research/brain_connectivity/data/precision/', 'test', n_iter);

load('output/matching_results/P_test_1_two_region_fixed_iter.mat')

% set convergence: epsilon = 0.01
% obj_{t} - obj_{t-1} < epsilon

n = 997;
conv_iter = zeros(n, 1);
eps = 0.001;
for i = 1:n
    % relative difference
    obj_diff = obj(i,2:n_iter) - obj(i,1:(n_iter-1)); % rdivide(obj(i,2:n_iter), obj(i,1:(n_iter-1))) - 1; 
    all_conv_iters = find(obj_diff < eps);
    conv_iter(i) = all_conv_iters(1);
end

conv_iter = conv_iter([1:(i_subj-1), (i_subj+1):n]);

tbl = tabulate(conv_iter);
x = categorical(tbl(:,1));
freq = tbl(:,3);
bar(x, freq)
xlabel('Number of iterations till convergence')
ylabel('Frequency among 996 mathces with subject 1')
title(sprintf('Stopping rule: obj(t)-obj(t-1) < eps=%.1e', eps))
saveas(gcf, 'figures/conv_eps_stopping_rule', 'jpg');


plot(0:(n_iter), obj(2:n, :))
title('Objective function: subj 1 - (subj 2-997)')
xlabel('Iterations')
ylabel('Values of objective function')
saveas(gcf, 'figures/conv_curve_all', 'jpg');

plot(0:(n_iter), obj(10, :))
title('Objective function: subj 1 - subj 10')
xlabel('Iterations')
ylabel('Values of objective function')
saveas(gcf, 'figures/conv_curve_1-10', 'jpg');



plot(0:(n_iter), obj(2:20, :))
title('Objective function: subj 1 - (subj 2-20)')
xlabel('Iterations')
ylabel('Values of objective function')
saveas(gcf, 'figures/conv_curve_1-2_20', 'jpg');


plot(0:(n_iter), offdiag_swap_counts_by_iter(2:20, :))
title('Number of swaps: subj 1 - (subj 2-20)')
xlabel('Iterations')
ylabel('Number of swaps')
saveas(gcf, 'figures/n_swaps_1-2_20', 'jpg');

plot(0:(n_iter), offdiag_swap_counts_by_iter(10, :))
title('Number of swaps: subj 1 - subj 10')
xlabel('Iterations')
ylabel('Values of objective function')
saveas(gcf, 'figures/n_swaps_1-10', 'jpg');

