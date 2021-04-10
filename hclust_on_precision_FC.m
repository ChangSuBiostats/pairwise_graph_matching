% Perform hierarchical clustering using precision FC maps
%
% Chang Su, c.su@yale.edu

% load precision FC maps
data_path = '/Users/chang/Documents/research/brain_connectivity/data/precision/';
addpath(data_path);

test = load('cc400_regFCprec_concat_hpf_997subj.mat');

% take median across 997 subjects
R = zeros(392);
subj_median = squeeze(median(test.C, 1));
% evaluate distance by F-norm
for i = 1:392
    for j = 1:392
        diff = subj_median(:, i) - subj_median(:, j); % subj1(:, i) - subj1(:, j);
        R(i, j) = dot(diff, diff);
    end
end
dist_mat = R;

% load ROI labels for left and right hemishperes
LR_labels = table2array(readtable('output/LR_labels.csv'));
left_index = find([LR_labels{:,2}]=='L');
right_index = find([LR_labels{:,2}]=='R');

% perform hierarchical clustering
dist_mat_right = dist_mat(right_index,right_index);
Z_right = linkage(dist_mat_right, 'complete');
dendrogram(Z_right);
title('complete linkage clustering, right hemishpere')
%saveas(gcf, 'figures/dendrogram_func_right.png')

dist_mat_left = dist_mat(left_index,left_index);
Z_left = linkage(dist_mat_left, 'complete');
dendrogram(Z_left);
title('complete linkage clustering, left hemishpere')
%saveas(gcf, 'figures/dendrogram_func_left.png')

% cut to have 11 clusters
n_cluster = 11;
T_left = cluster(Z_left,'maxclust',n_cluster);
T_right = cluster(Z_right,'maxclust',n_cluster);

% number of observations per cluster
tabulate(T_left)
tabulate(T_right)

% cluster size
left_count = histc(T_left, unique(T_left));
right_count = histc(T_right, unique(T_right));
cluster_count = reshape([left_count right_count], [2*n_cluster 1]);

% create an index vector for rearranging matrix
re_index_left = [left_index(find(T_left==1)) left_index(find(T_left==2)) left_index(find(T_left==3)) left_index(find(T_left==4)) ...
    left_index(find(T_left==5)) left_index(find(T_left==6)) left_index(find(T_left==7)) left_index(find(T_left==8)) ...
    left_index(find(T_left==9)) left_index(find(T_left==10)) left_index(find(T_left==11)) ];
re_index_right = [right_index(find(T_right==1)) right_index(find(T_right==2)) right_index(find(T_right==3)) right_index(find(T_right==4)) ...
    right_index(find(T_right==5)) right_index(find(T_right==6)) right_index(find(T_right==7)) right_index(find(T_right==8)) ...
    right_index(find(T_right==9)) right_index(find(T_right==10)) right_index(find(T_right==11)) ];
re_index = [re_index_left re_index_right];

% save the results
save('output/hclust_11_functional.mat','re_index', 'cluster_count', 'left_index', 'T_left', 'right_index', 'T_right');
