function [obj_err, Xs_rec, recovery_errs] = solver_gd(y, Bs, Ls_init, Rs_init, X_gt, s, n, r, K, maxiter, step_size, is_show)

if mod(n,2) == 0
    n1 = n/2;
    DD = [1:n1 n1 n1-1:-1:1].';
else
    n1 = (n+1)/2;
    DD = [1:n1 n1-1:-1:1].';
end
D = sqrt(DD);


obj_err = zeros(maxiter,1);
Ls = Ls_init;
Rs = Rs_init;
sn1 = size(Ls,1);
n2 = size(Rs,1);

recovery_errs = zeros(maxiter,1);

Zs_gt = zeros(sn1, n2, K);
for kk = 1:K
    Zs_gt(:,:,kk) = G(X_gt(:,:,kk)*diag(D));
end


s0 = zeros(K,1);
for kkk = 1:K
    s0(kkk) = norm(Ls(:,:,kkk));
end
s0 = max(s0(:));
% s0 = max(norm(Z0v),norm(Z0u));



for iter = 1:maxiter
    
    % choose step size via line search
    [curr_obj, grad_L, grad_R, grad_norm, ~, ~, ~] = get_obj_grad(y, Bs, Ls, Rs, s, n, r, K);
    obj_err(iter) = curr_obj;
    if ~step_size
        beta = s0^2;
        alpha = 0.5;
        [new_obj, ~, ~, ~, ~, ~, ~] = get_obj_grad(y, Bs, Ls, Rs, s, n, r, K);
        jj = 0;

        while new_obj > curr_obj - 1/100*beta /(s0^2) * grad_norm^2 || jj<=10
            jj = jj+1;
            beta = alpha * beta;
            Ls_tmp = Ls - beta/s0^2 * grad_L;
            Rs_tmp = Rs - beta/s0^2 * grad_R;
            [new_obj, ~, ~, ~,~,~,~] = get_obj_grad(y, Bs, Ls_tmp, Rs_tmp, s, n, r, K);
        end

    else
        beta = step_size;%/norm(Zs_gt(:), 'fro');
    end
    %fprintf('jj=%d\teta=%f\n',jj,beta);
    Ls = Ls - beta/s0^2 * grad_L;
    Rs = Rs - beta/s0^2 * grad_R;
    
    
    % 计算恢复误差
    Xs_rec = zeros(s,n,K);
    for kk = 1:K
        Zs = Ls(:,:,kk) * Rs(:,:,kk)';
        Xs_rec(:,:,kk) = Gstar(Zs,D, s)*diag(1./D);
        
    end
    recovery_errs(iter,1) = norm(Xs_rec(:) - X_gt(:), 'fro')/norm(X_gt(:), 'fro');
    if recovery_errs(iter,1)>=10^2
        recovery_errs = recovery_errs(1:iter,1);
        break;
    end
    if is_show
        fprintf('Iter=%d\tStep Size=%f\tLogObj=%f\tLogRela=%f\n', iter, beta, log10(obj_err(iter)), log10(recovery_errs(iter,1) ));
    end
end


end
