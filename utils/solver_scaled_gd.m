function [obj_err, Xs_rec, recovery_errs] = solver_scaled_gd(y, Bs, Ls_init, Rs_init, X_gt, s, n, r, K, maxiter, tol_rec,tol_gm, tol_obj, step_size, is_show)

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


Zs_gt = zeros(sn1, n2, K);
for kk = 1:K
    Zs_gt(:,:,kk) = G(X_gt(:,:,kk)*diag(D));
end


recovery_errs = zeros(maxiter,1);

for iter = 1:maxiter
    % choose step size via line search
    [curr_obj, ~, ~, ~, grad_s_L, grad_s_R, grad_s_norm] = get_obj_grad(y, Bs, Ls, Rs, s, n, r, K);
    obj_err(iter) = curr_obj;
    
    if ~step_size % 线搜索获得步长
        jj= 0;
        beta = 1;
        alpha = 0.5;
        Ls_tmp = Ls - beta * grad_s_L;
        Rs_tmp = Rs - beta * grad_s_R;
        [new_obj, ~, ~, ~, ~, ~, ~] = get_obj_grad(y, Bs, Ls_tmp, Rs_tmp, s, n, r, K);

        while new_obj > curr_obj - beta/2  * grad_s_norm^2 || jj<=5
            jj = jj+1;
            beta = alpha * beta;
            Ls_tmp = Ls - beta * grad_s_L;
            Rs_tmp = Rs - beta * grad_s_R;
            [new_obj, ~, ~, ~,~,~,~] = get_obj_grad(y, Bs, Ls_tmp, Rs_tmp, s, n, r, K);
        end
    else % 固定步长
        beta = step_size;
    end

    Ls = Ls - beta * grad_s_L;
    Rs = Rs - beta * grad_s_R;
    
    
    % 计算恢复误差
    Xs_rec = zeros(s,n,K);
    for kk = 1:K
        Zs = Ls(:,:,kk) * Rs(:,:,kk)';
        Xs_rec(:,:,kk) = Gstar(Zs,D, s)*diag(1./D);
        
    end
    recovery_errs(iter,1) = norm(Xs_rec(:) - X_gt(:), 'fro')/norm(X_gt(:), 'fro');
    %{
    %if recovery_errs(iter,1)>=10^2
        recovery_errs = recovery_errs(1:iter,1);
        break;
    end
    %}
    if is_show
        fprintf('Iter=%d\tStep Size=%f\tLogObj=%f\tLogRela=%f\n', iter, beta, log10(obj_err(iter)), log10(recovery_errs(iter,1)));
    end
    %{
    if X_err(iter) < tol_rec || grad_s_norm < tol_gm || curr_obj < tol_obj
        
        X_err = X_err(1:iter);
        for kk = 1:K
            Xs_rec(:,:,kk) = Gstar(Zs(:,:,kk), D, s)*diag(1./D);
        end

        return; 
    end
    %}

end






end
