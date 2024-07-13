clc; clear all; close all
%{
信号频率不可分

%}
n = 48;
K=2;

if mod(n,2) == 0
    n1 = n/2;
    DD = [1:n1 n1 n1-1:-1:1].';
else
    n1 = (n+1)/2;
    DD = [1:n1 n1-1:-1:1].';
end
n2 = n+1 - n1;
D = sqrt(DD);




empir_probs_scaledGD_r_vs_s_without_sep = zeros(6,6);

fileName = sprintf('data_empir_probs_scaledGD_r_vs_s_%d_%s_wihtout_sep.mat', n, datetime);


for  ii = 1:6
    for jj = 1:6
        r = ii; s = jj;
        sucind = 0;
        for avg = 1:20
            is_sep = 0;
            [X_gt, Bs, y] = get_data_demixing(n,s,K,r, is_sep);
            
            
            Zs_gt = zeros(s*n1, n2);
            for kk = 1:K
                Zs_gt(:,:,kk) = G(X_gt(:,:,kk)*diag(D));
            end
            maxiter = 200;
            tol_rec = 1e-5;
            tol_gm  = 1e-9;
            tol_obj = 1e-5;

            tic;
            [Ls_init, Rs_init, ~] = spectral_initial(y, Bs, n, s, r, K);
            [obj_scaled_gd_err, X_scaled_gd_err, Xs_rec] = solver_scaled_gd(y, Bs, Ls_init, Rs_init, Zs_gt, s, n, r, K, maxiter,tol_rec,tol_gm, tol_obj);
            time_cost = toc;
            relerror =  norm(Xs_rec(:) - X_gt(:))/norm(X_gt(:));
            fprintf('r=%d\ts=%d\tErr=%.4f\tTime=%.4f\n', r, s, relerror, time_cost);
            
            %Data_name = sprintf('phaseTransition_scaled_GD_%d.txt', n);
            %fid = fopen(Data_name,'a');
            %fprintf(fid,'r=%d\ts=%d\tErr=%.4f\n', r, s, relerror);
            %fclose(fid);
            
            
            if relerror <= 1e-3
                sucind = sucind+1;
            end
        end
        empir_probs_scaledGD_r_vs_s_without_sep(ii,jj) = sucind/20; 
        if empir_probs_scaledGD_r_vs_s_without_sep(ii,jj) ==0
            for jjj = jj+1:6
                empir_probs_scaledGD_r_vs_s_without_sep(ii,jjj) = 0;
            end
            break;
        end
    end
end

save(fileName, 'empir_probs_scaledGD_r_vs_s_without_sep')
