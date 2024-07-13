function [fs_rec,ps, peaks] = MUSIC(X, r, grid_size)
% this function is the realization of music algorith
% Input:
% X: is the data matrix, the column space is spaned by {a_fk}
% [n, s] = size(X)
% where n is the dim and s is the number of snapshots
% r is the number of spikes
% grid_size: 

[n,~] = size(X);
[U, ~, ~] = svd(X);
En=U(:,r+1:end);    

f = 0:grid_size:1;
ps = zeros(length(f),1);
for kk = 1:length(f)
    ps(kk) = pseudospectrum(f(kk), En, n);
end
[~, fs_res_loc] = findpeaks(ps, 'NPeaks', r, 'SortStr', 'descend');

ts = 0:grid_size:1;

fs_rec = ts(fs_res_loc);
peaks = ps(fs_res_loc);

end