clc; clear all; close all
%{
simultaneous blind super-resolution and demixing
这个代码说明scaled-gradient descent 线性收敛，并且收敛速率不依赖于条件数

%}
n = 128;
r = 2;
s = 2;
K = 2;

if mod(n,2) == 0
    n1 = n/2;
    DD = [1:n1 n1 n1-1:-1:1].';
else
    n1 = (n+1)/2;
    DD = [1:n1 n1-1:-1:1].';
end
n2 = n+1 - n1;
D = sqrt(DD);



%% 超参数
max_iter = 1000;
tol_rec = 1e-4;
tol_gm  = 1e-8;
tol_obj = 1e-5;
kappas = [1, 5 10 20];

scaled_gd_recover_errs = zeros(max_iter, length(kappas));
gd_recover_errs = zeros(max_iter, length(kappas));
for ka = 1:length(kappas)
    X_gt = zeros(s,n,K);
    Bs = zeros(n,s,K);
    y = zeros(n,1);

    kappa = kappas(ka);
    for kk = 1:K
        [~, ~, ~, ~, X_gt(:,:,kk)] = get_X_with_sep_fixed_condition_number(r, s, n, kappa);
        Bs(:,:,kk) = -sqrt(3)+2*sqrt(3).*rand(n,s);
        y = y + diag(Bs(:,:,kk)*X_gt(:,:,kk));
    end




    %% 初始值
    [Ls_init, Rs_init, Xs_init] = spectral_initial(y, Bs, n,s, r, K);
    %% scaled gradient descent

    step_size = 0.6/K; %1 常数步长; 0 线搜索步长；
    [~,~, scaled_gd_recover_errs(:,ka)] = solver_scaled_gd(y, Bs, Ls_init, Rs_init, X_gt, s, n, r, K, max_iter, tol_rec,tol_gm, tol_obj, step_size);
    
    %% gradient descent
    step_size = 0.6/K;
    [~,~, gd_recover_errs(:,ka)] = solver_gd(y, Bs, Ls_init, Rs_init, X_gt, s, n, r, K, max_iter, step_size);
end

maker_idx = 1:10:max_iter;
semilogy(1:max_iter, scaled_gd_recover_errs(:,1), '-o', 'MarkerIndices',maker_idx, 'LineWidth', 2.0);
%plot(1:max_iter, log10(scaled_gd_recover_errs(:,1)), '-o', 'MarkerIndices',maker_idx, 'LineWidth', 2.0);

hold on;
semilogy(1:max_iter, scaled_gd_recover_errs(:,2), '-o','MarkerIndices',maker_idx, 'LineWidth', 2.0);


semilogy(1:max_iter, scaled_gd_recover_errs(:,3), '-o','MarkerIndices',maker_idx, 'LineWidth', 2.0);
semilogy(1:max_iter, scaled_gd_recover_errs(:,4), '-o','MarkerIndices',maker_idx, 'LineWidth', 2.0);
%legend('sgd 1', 'sgd 5', 'sgd 10', 'sgd 20');
%}

semilogy(1:max_iter, gd_recover_errs(:,1), '-^','MarkerIndices',maker_idx, 'LineWidth', 2.0);

semilogy(1:max_iter, gd_recover_errs(:,2), '-^','MarkerIndices',maker_idx, 'LineWidth', 2.0);

semilogy(1:max_iter, gd_recover_errs(:,3), '-^','MarkerIndices',maker_idx, 'LineWidth', 2.0);
semilogy(1:max_iter, gd_recover_errs(:,4), '-^','MarkerIndices',maker_idx, 'LineWidth', 2.0);

legend('sgd 1', 'sgd 5', 'sgd 10', 'sgd 20', 'gd 1', 'gd 5', 'gd 10', 'gd 20');
%}


