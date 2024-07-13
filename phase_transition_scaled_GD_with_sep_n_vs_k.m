clc; clear all; close all

s = 3;
r = 3;



%ns = 80:10:160;
ns = 170:10:240;
Ks = 2:1:8;
empir_probs_scaledGD_n_vs_K_2 = zeros(6,6);

fileName = sprintf('data_empir_probs_scaledGD_n_vs_K_s%d_r%d_%s.mat', s,r, datetime);


for  ii = 1:length(ns)
    for jj = 1:length(Ks)
        n = ns(ii);
        K = Ks(jj);

        if mod(n,2) == 0
            n1 = n/2;
            DD = [1:n1 n1 n1-1:-1:1].';
        else
            n1 = (n+1)/2;
            DD = [1:n1 n1-1:-1:1].';
        end
        n2 = n+1 - n1;
        D = sqrt(DD);


        sucind = 0;
        for avg = 1:20
            is_sep = 1;
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
            fprintf('n=%d\tK=%d\tErr=%.4f\tTime=%.4f\n', n, K, relerror, time_cost);
            
            %Data_name = sprintf('phaseTransition_scaled_GD_%d.txt', n);
            %fid = fopen(Data_name,'a');
            %fprintf(fid,'r=%d\ts=%d\tErr=%.4f\n', r, s, relerror);
            %fclose(fid);
            
            
            if relerror <= 1e-3
                sucind = sucind+1;
            end
        end
        empir_probs_scaledGD_n_vs_K_2(ii,jj) = sucind/20; 
        if empir_probs_scaledGD_n_vs_K_2(ii,jj) ==0
            for jjj = jj+1:length(Ks)
                empir_probs_scaledGD_n_vs_K_2(ii,jjj) = 0;
            end
            break;
        end
    end
end

save(fileName, 'empir_probs_scaledGD_n_vs_K_2')
