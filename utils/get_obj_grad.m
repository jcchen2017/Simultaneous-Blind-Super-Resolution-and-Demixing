function [obj, grad_L, grad_R, grad_norm, grad_s_L, grad_s_R, grad_s_norm] = get_obj_grad(y, Bs, Ls, Rs, s, n, r, K)
%{
给定一点，计算该点的目标函数值以及梯度

%}

if mod(n,2) == 0
    n1 = n/2;
    DD = [1:n1 n1 n1-1:-1:1].';
else
    n1 = (n+1)/2;
    DD = [1:n1 n1-1:-1:1].';
end
n2 = n+1 - n1;
weight_d = sqrt(DD);


Zs = zeros(s*n1, n2, K);
for kk = 1:K
    Zs(:,:,kk) = Ls(:,:,kk) * Rs(:,:,kk)';
end


res1 = zeros(n,1);

for kk = 1:K
    GstarZ = Gstar(Zs(:,:,kk), weight_d, s); % s * n; 
    res1 = res1 + diag(Bs(:,:,kk)*GstarZ);
end
res1 = res1 - diag(weight_d) * y; %n*1

obj1 = 0.5*norm(res1)^2;
obj2 = 0;
proj_Z = zeros(s*n1, n2, K);

% computer gradients
grad1_L = zeros(s*n1, r, K);
grad1_R = zeros(n2, r, K);
grad2_L = zeros(s*n1, r, K);
grad2_R = zeros(n2, r, K);


for kk = 1:K
    GstarZ = Gstar(Zs(:,:,kk), weight_d, s);
    proj_Z(:,:,kk) = Zs(:,:,kk) - G(GstarZ);
    obj2 = obj2 + 0.5*norm(proj_Z(:,:,kk), 'fro')^2;
    
    tmp_grad1 = G(Bs(:,:,kk)' * diag(res1));
    grad1_L(:,:,kk) = tmp_grad1 * Rs(:,:,kk);
    grad1_R(:,:,kk) = tmp_grad1' * Ls(:,:,kk);
    grad2_L(:,:,kk) = proj_Z(:,:,kk) * Rs(:,:,kk);
    grad2_R(:,:,kk) = (proj_Z(:,:,kk))' * Ls(:,:,kk);
    
    
end

obj = obj1 + obj2;
grad_L = grad1_L + grad2_L;
grad_R = grad1_R + grad2_R;

% 计算梯度的范数

grad_norm = 0;
for kk = 1:K
    gm_k = norm(grad_L(:,:,kk), 'fro')^2 + norm(grad_R(:,:,kk), 'fro')^2;
    grad_norm = grad_norm + gm_k;
end
grad_norm = sqrt(grad_norm);

grad_s_L = zeros(s*n1, r, K);
grad_s_R = zeros(n2, r, K);
% 计算scaled gradient
for kk = 1:K
    RtR = Rs(:,:,kk)' * Rs(:,:,kk);
    LtL = Ls(:,:,kk)' * Ls(:,:,kk);
    grad_s_L(:,:,kk) = grad_L(:,:,kk) * inv(RtR);
    grad_s_R(:,:,kk) = grad_R(:,:,kk) * inv(LtL);
end

grad_s_norm = 0;
for kk = 1:K
    gm_s_k = norm(grad_s_L(:,:,kk), 'fro')^2 + norm(grad_s_R(:,:,kk), 'fro')^2;
    grad_s_norm = grad_s_norm + gm_s_k;
end
grad_s_norm = sqrt(grad_s_norm);


end