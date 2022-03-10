function [P_old, diff, obj, new_m, cost_m, sum_swaps] = iterative_procedure(m1, m2, n_iter, lambda, penalty_m)
% match two maps using the iterative procedure, with possile regularization by penalizing certain matches
% Arguments:
%   m1: map 1
%   m2: map 2, to be matched with map 1
%   n_iter: number of iterations in the iterative procedure
%   lambda: penalty weight
%   penalty_m: penalty matrix specifying the penalties for matches between ROIs
% Return:
%   P_old: permutation matrix
%   diff: a vector of Forbenius norm between original and permuted maps,
%   along iterations
%   obj: values of objective function along iterations, = diff + lambda * trace(P'D)
%   new_m: permuted m1, to match with m2

    n_dim = size(m1,1);

    %% normalize two maps before matching
    m1 = m1 / norm(m1, 'fro');
    m2 = m2 / norm(m2, 'fro');
    
    diff = zeros(1,n_iter+1);
    obj = zeros(1,n_iter+1);
    sum_swaps = zeros(1, n_iter+1);
                                                                                
    %% initiate the iterative procedure
    
    P_old = eye(n_dim);

    %% objective function
    diff(1,1) = norm(P_old' * m1 * P_old - m2,'fro')^2;
    obj(1,1) = trace(P_old' * (m1' * P_old * m2 - lambda * penalty_m));
    sum_swaps(1,1) = 0;
    
    %% iterative updates
    % tic
    % stoppint criterion
    eps = 0.001;
    rdelta = 1;
    iter = 1;
    % for iter = 1:n_iter
    while rdelta > eps && iter <= n_iter
        % -------------------
        M = matchpairs(m1' * P_old * m2 - lambda * penalty_m, -99999, 'max');
        P_new = full(sparse(M(:, 1), M(:, 2), 1, n_dim, n_dim)); % adjust the format
        obj(1, iter + 1) = trace(P_new' * (m1' * P_old * m2 - lambda * penalty_m));
        P_old = P_new;
        % -------------------
        diff(1, iter + 1) = norm(P_old' * m1 * P_old - m2, 'fro')^2;
        tmp = P_old - diag(diag(P_old));
        sum_swaps(1, iter + 1) = sum(tmp, 'all');
        % stopping rule: difference <= epsilon
        rdelta = obj(1, iter + 1) - obj(1, iter);
        iter = iter + 1;
    end
    % toc
    
    if iter > n_iter
        fprintf('Iterative proceure did not converge in %i iterations.', i_ref)    
    end

    new_m = P_old' * m1 * P_old;
    
    % new: evaluate the cost matrix used in the Hungarian algorithm
    cost_m = m1' * P_old * m2 - lambda * penalty_m;
