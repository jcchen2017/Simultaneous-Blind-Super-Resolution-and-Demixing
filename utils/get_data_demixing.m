function [X_gt, Bs, y] = get_data_demixing(n,s,K,r, is_sep)

X_gt = zeros(s,n,K);
Bs = zeros(n,s,K);
y = zeros(n,1);
dynamic_range = 20;
for kk = 1:K
    if is_sep
    [~, ~, ~, ~, X_gt(:,:,kk)] = get_X_with_sep(r, s, n, dynamic_range);
    else
    [~, ~, ~, ~, X_gt(:,:,kk)] = get_X_without_sep(r, s, n, dynamic_range);
    end
    %Bs(:,:,kk) = get_random_FFT_B(n,s);
    %Bs(:,:,kk) = randn(n,s);
    %Bs(:,:,kk) = 2*binornd(1,0.5, n,s)-1;
    Bs(:,:,kk) = -sqrt(3)+2*sqrt(3).*rand(n,s);
    y = y + diag(Bs(:,:,kk)*X_gt(:,:,kk));
end


end
