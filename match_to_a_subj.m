function []=match_to_a_subj(i_ref, lambda, penalty, data_path)
    arguments
        i_ref (1,1) {mustBeNumeric,mustBeReal} = 1
        lambda (1,1) {mustBeNumeric,mustBeReal} = 3e-4
        penalty (1,:) char {mustBeMember(penalty,{'vanilla','two_region'})} = 'two_region'
        data_path (1,:) char = 'data/'
    end
    
    fprintf('i_ref=%i, lambda=%.1e, penalty=%s, data_path=%s \n', i_ref, lambda, penalty, data_path)
    
    % reference: 
    % 1. make MATLAB take arguments from command line:
    %   https://stackoverflow.com/questions/8981168/running-a-matlab-program-with-arguments
    % 2. set default value for parameters
    %   https://www.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html
    % 3. save large sparse matrix efficiently
    %   reference: https://www.mathworks.com/matlabcentral/answers/101798-how-can-i-write-a-sparse-matrix-s-non-zero-values-and-the-corresponding-row-and-column-information-t
        
    
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
    % retest = load('cc400_regFCprec_concat_hpf_retest41.mat');
    n = size(test.subj, 2);
    
    %% start matching
    
    % set penalty matrices
    if strcmp(penalty, 'vanilla')
        penalty_m = load('output/vanilla_penalty.mat').vanilla_p;
    elseif strcmp(penalty, 'two_region')
        penalty_m = load('output/two_regions_penalty.mat').two_region_p;
    end    
    
    % set the FC map from subj i as the reference map
    m_ref = squeeze(test.C(i_ref, :, :));
    m_ref(logical(eye(392))) = 0; % remove the effect of diagonal values
    
    % set number of iterations to 1 to save computational time
    n_iter = 1;
    
    % keep track of
    % 1) the position of nonzero entries in the permutation matrix
    swap_positions = zeros(n, 392, 2);
    % 2) difference between two matrices before and after permutation
    diff = zeros(n,2); 
    % 3) number of off diagonal swaps
    offdiag_swap_counts = zeros(n,1); 
    % 4) array of permutation matrices 
    P_array = zeros(392, 392, n);
    
    % Loop through all pairs 
    % Though the computation is redundant (some pairs might have been
    % matched before),
    % we stick with this because the data will be eaiser to access.
    for i = 1:n
        m2 = squeeze(test.C(i, :, :));
        m2(logical(eye(392))) = 0;
        [P, d, ~, ~] = iterative_procedure(m_ref, m2, n_iter, lambda, penalty_m);
        diff(i,:) = d;
        offdiag_swap_counts(i) = 392 - sum(P(logical(eye(392))));
        [row, col, ] = find(P);
        % temporarily store the permutation matrices
        P_array(:, :, i) = P;
        swap_positions(i, :, :) = [row, col];
    end
    
    % sum of permutation matrices
    [row, col, v] = find(sum(P_array, 3));
    % the positions and values of nonzero entries
    sum_swaps = [row, col, v];

    save(strcat('output/matching_results/P_', num2str(i_ref), '.mat'), 'swap_positions', 'sum_swaps', 'diff', 'offdiag_swap_counts');
end

