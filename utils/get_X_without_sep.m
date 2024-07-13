function [fs, cs, H, A, X0] = get_X_without_sep(r, s, n, dynamic_range)

% frequency
fs = rand(r,1);
%{
fs_sort = sort(fs);
if r==1
    sep = 1;
else
    sep = min(min(abs(fs_sort(2:end)-fs_sort(1:end-1))),fs_sort(1)+1-fs_sort(end));
end
while (sep<1/n)
    fs = rand(r,1);
    fs_sort = sort(fs);
    if r==1
        sep = 1;
    else
        sep = min(min(abs(fs_sort(2:end)-fs_sort(1:end-1))),fs_sort(1)+1-fs_sort(end));
    end
end
%}

% amplitude
dynamic_range= 20;
cs = exp(-1i*2*pi*rand(r,1)) .* (1 + 10.^(rand(r,1).*(dynamic_range/20)));
%cs = exp(-1i*2*pi*rand(r,1)) .* linspace(1/dynamic_range,1, r)';

% ill-conditioned 
%cs = cs.*linspace(1/10,1,r)';

Cs = diag(cs); %J x J diag matrix

H = zeros(s, r);
for kk = 1:r
    H(:,kk) = complex(randn(s,1), randn(s,1))/sqrt(2);
    H(:,kk) = H(:, kk)/norm(H(:,kk));
end

A = zeros(r, n);
for kk  = 1:r
    tmp = exp(1i*2*pi*fs(kk).*(0:1:n-1))/sqrt(n); % r* n
    A(kk,:) = tmp';
end
X0 = H * Cs * A;


end



