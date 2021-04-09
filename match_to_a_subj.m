function []=match_to_a_subj(i_ref, lambda, penalty)
    arguments
        i_ref (1,1) {mustBeNumeric,mustBeReal}
        lambda (1,1) {mustBeNumeric,mustBeReal} = 3e-4
        penalty (1,1) char {mustBeMember(penalty,{'vanilla','two_region'})} = 'vanilla'
    end
    
    % reference: 
    % make MATLAB take arguments from command line:
    %   https://stackoverflow.com/questions/8981168/running-a-matlab-program-with-arguments
    % set default value for parameters
    %   https://www.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html
    
    % load codes for graph matching
    addpath /Users/chang/Documents/research/brain_connectivity/code
    
    % load data
    test = load('/Users/chang/Documents/research/brain_connectivity/data/precision/cc400_regFCprec_concat_hpf_997subj.mat');
    retest = load('/Users/chang/Documents/research/brain_connectivity/data/precision/test_retest/cc400_regFCprec_concat_hpf_retest41.mat');
    n = size(test.subj, 2);
    
    % design a penalty matrix
    if penalty == 'vanilla'
    
    % set the FC map from subj i as the reference map
    m_ref = squeeze(test.C(i_ref, :, :));
    m_ref(logical(eye(392))) = 0;
    
    for i = 1:n
        m2 = squeeze(test.C(i, :, :));
        m2(logical(eye(392))) = 0;
        [P, diff, obj, new_m] = iterative_procedure(m_ref, m2, n_iter, lambda, vanilla_p);
        P_cross{1,i} = P;
        diff_m_cross(i) = diff(n_iter + 1);
        off_diag_matches_cross(1, i) = 392 - sum(P_cross{1, i}(logical(eye(392))));
    end


end