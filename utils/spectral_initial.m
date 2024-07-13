function [Ls_init, Rs_init, Xs_init] = spectral_initial(y, Bs, n,s, r, K)
%{
通过谱方法获得初始值
y: n*1 观测结果
A: s*n -> n 采样算子

%}
if mod(n,2) == 0
    n1 = n/2;
    DD = [1:n1 n1 n1-1:-1:1].';
else
    n1 = (n+1)/2;
    DD = [1:n1 n1-1:-1:1].';
end
n2 = n+1 - n1;
D = sqrt(DD);
Ls_init = zeros(s*n1, r, K);
Rs_init = zeros(n2, r, K);
Xs_init = zeros(s, n, K);
for kk = 1:K
    
    Aty = Bs(:,:,kk)'*diag(y); % s*n
    HAty = zeros(s*n1, n2);
    for j1 = 1:n1
        for j2 = 1:n2
            row_idx = (j1-1)*s+1:j1*s;
            HAty(row_idx, j2) = Aty(:, j1+j2-1);
        end
    end

    [U0, S0, V0] = svds(HAty, r);
    Ls_init(:,:,kk) = U0*sqrt(S0);
    Rs_init(:,:,kk) = V0*sqrt(S0);
    Xs_init(:,:,kk) = Gstar(Ls_init(:,:,kk) * Rs_init(:,:,kk)',D, s)*diag(1./D);
end



end