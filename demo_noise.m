clc; clear all; close all
%{
simultaneous blind super-resolution and demixing
噪声下恢复性能测试

%}
% 
% n =256;r=s=2;K=4;
n = 512; 
r = 4;
s = 4;
K = 6;

if mod(n,2) == 0
    n1 = n/2;
    DD = [1:n1 n1 n1-1:-1:1].';
else
    n1 = (n+1)/2;
    DD = [1:n1 n1-1:-1:1].';
end
n2 = n+1 - n1;
D = sqrt(DD);


avgTimes = 20;

snrs = 0:10:60;
kappa = 15; %1,5, 15
recovery_err_noise = zeros(length(snrs),avgTimes);
for idx_snr = 1:length(snrs)
    parfor mc = 1:avgTimes
        
        snr = snrs(idx_snr);
        fprintf('SNR=%fdB\tAvg=%d\n', snr,mc);
        % get data
        X_gt = zeros(s,n,K);
        Bs = zeros(n,s,K);
        y = zeros(n,1);

        kappa = 1;
        for kk = 1:K
            [~, ~, ~, ~, X_gt(:,:,kk)] = get_X_with_sep_fixed_condition_number(r, s, n, kappa);
            Bs(:,:,kk) = -sqrt(3)+2*sqrt(3).*rand(n,s);
            y = y + diag(Bs(:,:,kk)*X_gt(:,:,kk));
        end
        ym = norm(y(:));
        % additive noise
        sigma = ym*10^(-snr/20);
        z = complex(randn(n,1), randn(n,1))/sqrt(2);
        z = sigma * z/norm(z(:));
        y = y+z;

        %% ground truth used to computer recovery error
        Ls_gt = zeros(s*n1, r,K);
        Rs_gt = zeros(n2, r, K);
        Zs_gt = zeros(s*n1, n2);
        for kk = 1:K
            Zs_gt(:,:,kk) = G(X_gt(:,:,kk)*diag(D));
            [Ls, Ss, Rs] = svds(Zs_gt(:,:,kk), r);
            Ls_gt(:,:,kk) = Ls* sqrt(Ss);
            Rs_gt(:,:,kk) = Rs * sqrt(Ss);
        end

        % 初始值
        [Ls_init, Rs_init, Xs_init] = spectral_initial(y, Bs, n,s, r, K);


        %% gradient descent
        maxiter = 400;
        tol_rec = 1e-4;
        tol_gm  = 1e-8;
        tol_obj = 1e-5;


        is_show = 0;

        step_size = 0.6/K; %1 常数步长; 0 线搜索步长；
        [~,~, scaled_gd_recover_err] = solver_scaled_gd(y, Bs, Ls_init, Rs_init, X_gt, s, n, r, K, maxiter, tol_rec,tol_gm, tol_obj, step_size, is_show);

        recovery_err_noise(idx_snr, mc) = scaled_gd_recover_err(end);
    end
end



plot(snrs, log10(mean(recovery_err_noise,2)),'-o');
legend('scaled gd');
%figure()
%plot(1:maxiter, log10(X_scaled_gd_err));