clc; clear all; close all

n = 48;
K=2;
empir_probs_vhl_r_vs_s = zeros(6,6);

fileName = sprintf('data_empir_probs_vhl_r_vs_s_%d_%s.mat', n, datetime);


for  i = 1:6
    for j = 1:6
        r = i; s = j;
        sucind = 0;
        for iter = 1:20
            %[~, ~ , ~,~,  ~,~,  ~,~,  X01, X02, B, y,D] = getSignals_wang3(r, r, s, n);
            [X_gt, Bs, y] = get_data_demixing(n,s,K,r);
            
            %[X1 , X2,  ~, ~ ] = solvermultiHankel(y, B, D);
            [X1 , X2] = solver_vhl(y, Bs, s, n);
            
            relerror =  sqrt((norm(X1 - X_gt(:,:,1), 'fro')^2+norm(X2-X_gt(:,:,2), 'fro')^2)/norm(X_gt(:))^2);
            fprintf('r=%d\ts=%d\tErr=%.4f\n', r, s, relerror);
            
            Data_name = sprintf('phaseTransition_%d.txt', n);
            fid = fopen(Data_name,'a');
            fprintf(fid,'r=%d\ts=%d\tErr=%.4f\n', r, s, relerror);
            fclose(fid);
            
            
            if relerror < 1e-6
                sucind = sucind+1;
            end
        end
        empir_probs_vhl_r_vs_s(i,j) = sucind/20; 
        if empir_probs_vhl_r_vs_s(i,j) <=1e-4
            for jj = j+1:6
                empir_probs_vhl_r_vs_s(i,j) = 0;
            end
            break;
        end
    end
end

save(fileName, 'empir_probs_vhl_r_vs_s')
colormap(gray);
imagesc(empir_probs_vhl_r_vs_s);
