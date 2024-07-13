function [fs1, fs2, cs1,cs2, H1,H2, A1,A2, X01, X02, B, y] = getSignals_with_sep(r1,r2, s, n)

fs1 = rand(r1,1);
fs2 = rand(r2,1);

if r1 == 1
    sep = 1;
else
    fs_sort = sort(fs1);
    sep = min(min(abs(fs_sort(2:end)-fs_sort(1:end-1))),fs_sort(1)+1-fs_sort(end));
end

while(sep<1/n)
    fs1 = rand(r1,1);
    if r1 == 1
        sep = 1;
    else
        fs_sort = sort(fs1);
        sep = min(min(abs(fs_sort(2:end)-fs_sort(1:end-1))),fs_sort(1)+1-fs_sort(end));
    end
end

dynamic_range = 20;
cs1 = exp(-1i*2*pi*rand(r1,1)) .* (1 + 10.^(rand(r1,1).*(dynamic_range/20)));
cs2 = exp(-1i*2*pi*rand(r2,1)) .* (1 + 10.^(rand(r2,1).*(dynamic_range/20)));


Cs1 = diag(cs1); %J x J diag matrix
Cs2 = diag(cs2);

H1 = zeros(s, r1);
H2 = zeros(s,r2);
for kk = 1:r1
    H1(:,kk) = randn(s,1);
    H1(:,kk) = H1(:, kk)/norm(H1(:,kk));
end

for kk = 1:r2
    H2(:,kk) = randn(s,1);
    H2(:,kk) = H2(:, kk)/norm(H2(:,kk));
end

A1 = zeros(r1, n); A2 = zeros(r2,n);
for kk  = 1:r1
    tmp = amplitude(fs1(kk), n);
    A1(kk,:) = tmp';
end
for kk  = 1:r2
    tmp = amplitude(fs2(kk), n);
    A2(kk,:) = tmp';
end
X01 = H1 * Cs1 * A1;
X02 = H2 * Cs2 * A2;


% get the observation

% Fourier samples
F1 = fft(eye(s));
B{1} = zeros(n,s);
for i = 1 : n
    B{1}(i,:) = F1(randi(s),:);
end

F2 = fft(eye(s));
B{2} = zeros(n,s);
for i = 1 : n
    B{2}(i,:) = F2(randi(s),:);
end

y = zeros(n,1);

for i = 1:n
    y(i,1) =  B{1}(i,:) * X01(:,i) + B{2}(i,:) * X02(:,i);
end

end

%%
function a = amplitude(fs, n)
N = 0:1:n-1;
a = exp(1i*2*pi*fs.*N);
end