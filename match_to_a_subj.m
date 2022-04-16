function []=match_to_a_subj(i_ref, lambda, penalty, data_path, ref_popu, n_iter)
%DOCFUN match a selected subject to all 997 subjects 
%
% Perform matching with precision FC maps
% Save the permutation matrices and relevant summary statistics
%
% INPUT
%   i_ref: the selected subject, which will be matched to all 997 subjects.
%           Default to 1.
%   lambda: weight for the penalty. Default to 3e-4.
%   penalty: name of the penalty matrix. Currently we support 'vanilla'
%       (penalizing all swaps), 'two_region' (penalizing only interactions across
%       two hclust blocks), 'yeo' (penalizing blocks within and between
%       Limbic, subcortical and cerebellum), 'yeo_all_pairs' (penalizing
%       all pairs that have at least one region in the Limbic, subcortical
%       or cerebellum networks under the yeo mapping. Default to
%       'yeo_all_pairs'
%   data_path: the path to HCP precision matrices 
%       cc400_regFCprec_concat_hpf_997subj.mat
%       cc400_regFCprec_concat_hpf_retest41.mat
%   ref_popu: which population does the reference matrix come from, i.e.
%       test or re-test. Default to test.
%   n_iter: maximum number of iterations to run; stops when the stopping
%       rule is met
%
% OUTPUT
%   a matlab file at output/matching_results/, named 'P_%i.mat' % i_ref,
%   which contains 
%       swap_positions: (997*392*2), each row corresponds to the
%       coordinates of the swaps
%       sum_swaps: (392,3), the frequency of swaps when matching to 997
%           subjects, where the first two columns are coordinates and the
%           thrid is frequency counts.
%       diff: (997,2), squared F-norm of the difference between two maps before and
%           after matching
%       offdiag_swap_counts: (997,1), number of offdiag swaps
%       obj: (997, n_iter+1), objective function across iterations
%       offdiag_swap_counts_by_iter: (997, n_iter+1), number of swaps across iterations
%
% Examples:
%   matlab -nodisplay -nodesktop -r "match_to_a_subj(1,3e-4,'yeo', \
%   '/Users/chang/Documents/research/brain_connectivity/data/precision/', 'test', 1)"
%
% Author:      Chang Su
% email:        c.su@yale.edu
% Matlab ver.:  R2019b
% Date:         10-Apr-2021
%
% Documentation reference: https://www.mathworks.com/matlabcentral/fileexchange/41423-documenting-help-section-of-an-m-file
%
    arguments
        i_ref (1,1) {mustBeNumeric,mustBeReal} = 1
        lambda (1,1) {mustBeNumeric,mustBeReal} = 3e-4
        penalty (1,:) char {mustBeMember(penalty,{'vanilla','two_region','yeo', 'yeo_all_pairs'})} = 'yeo_all_pairs'
        data_path (1,:) char = 'data/'
        ref_popu (1,:) char = 'test'
        n_iter (1,1) {mustBeNumeric,mustBeReal} = 1
    end
    
    fprintf('i_ref=%i, lambda=%.1e, penalty=%s, data_path=%s , ref_popu=%s, n_iter=%i \n', i_ref, lambda, penalty, data_path, ref_popu, n_iter)    
    
    %% set up paths and directories
    % add the path where data are stored
    addpath(data_path)
    
    if ~exist('output', 'dir')
       fprintf('directory ./output created \n')
       mkdir('output')
    end
    
    if ~exist('output/matching_results', 'dir')
       fprintf('directory ./output/matching_results created \n')
       mkdir('output/matching_results')
    end
    
    %% load data
    test = load('cc400_regFCprec_concat_hpf_997subj.mat');
    if strcmp(ref_popu, 'retest')
        retest = load('cc400_regFCprec_concat_hpf_retest41.mat');
    end
    n = size(test.subj, 2);
    
    %% start matching
    
    % set penalty matrices
    if strcmp(penalty, 'vanilla')
        penalty_m = load('output/vanilla_penalty.mat').vanilla_p;
    elseif strcmp(penalty, 'two_region')
        penalty_m = load('output/two_regions_penalty.mat').two_region_p;
    elseif strcmp(penalty, 'yeo')
        penalty_m = load('output/yeo_penalty.mat').yeo_p;
    elseif strcmp(penalty, 'null')
        penalty_m = load('output/null_penalty.mat').null_p;
    elseif strcmp(penalty, 'yeo_all_pairs')
        penalty_m = load('output/yeo_all_pairs_penalty.mat').yeo_all_pairs;
    end    
    
    % set the FC map from subj i in ref_popu as the reference map
    if strcmp(ref_popu, 'test')
        m_ref = squeeze(test.C(i_ref, :, :));
    elseif strcmp(ref_popu, 'retest')
        m_ref = squeeze(retest.C(i_ref, :, :));
    end
    m_ref(logical(eye(392))) = 0; % remove the effect of diagonal values
    
    % set number of iterations to 1 to save computational cost
    % n_iter = 1;
    
    % keep track of
    % 1) the position of nonzero entries in the permutation matrix
    swap_positions = zeros(n, 392, 2);
    % 2) difference between two matrices before and after permutation
    diff = zeros(n,n_iter+1); 
    % 3) number of off diagonal swaps
    offdiag_swap_counts = zeros(n,1); 
    % 4) array of permutation matrices 
    P_array = zeros(392, 392, n);
    % 5) objective functions
    obj = zeros(n,n_iter+1);
    % 6) sum of swaps
    offdiag_swap_counts_by_iter = zeros(n, n_iter+1);
    
    % Loop through all pairs 
    % Though the computation is redundant (some pairs might have been
    % matched before),
    % we stick with this because the data will be eaiser to access.
    for i = 1:n
        m2 = squeeze(test.C(i, :, :));
        m2(logical(eye(392))) = 0;
        [P, d, ob, ~, ~, ss] = iterative_procedure(m_ref, m2, n_iter, lambda, penalty_m);
        diff(i,:) = d;
        obj(i,:) = ob;
        offdiag_swap_counts(i) = 392 - sum(P(logical(eye(392))));
        offdiag_swap_counts_by_iter(i,:) = ss;
        [row, col, ] = find(P);
        % temporarily store the permutation matrices
        P_array(:, :, i) = P;
        swap_positions(i, :, :) = [row, col];
    end
    
    % sum of permutation matrices
    [row, col, v] = find(sum(P_array, 3));
    % the positions and values of nonzero entries
    sum_swaps = [row, col, v];

    save(strcat('output/matching_results/P_', ref_popu, '_', num2str(i_ref), '_', penalty, '.mat'), 'swap_positions', 'sum_swaps', 'diff', 'offdiag_swap_counts', 'obj', 'offdiag_swap_counts_by_iter');
end

% reference: 
    % 1. make MATLAB take arguments from command line:
    %   https://stackoverflow.com/questions/8981168/running-a-matlab-program-with-arguments
    % 2. set default value for parameters
    %   https://www.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html
    % 3. save large sparse matrix efficiently
    %   reference: https://www.mathworks.com/matlabcentral/answers/101798-how-can-i-write-a-sparse-matrix-s-non-zero-values-and-the-corresponding-row-and-column-information-t
    
