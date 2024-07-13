clear all;close all;clc

addpath('./utils/');

n = 256; 
r = 4;
s = 4;
K = 2;
if mod(n,2) == 0
    n1 = n/2;
    D = [1:n1 n1 n1-1:-1:1].';
else
    n1 = (n+1)/2;
    D = [1:n1 n1-1:-1:1].';
end
n2 = n+1 - n1;

% get data
X_gt = zeros(s,n,K);
Bs = zeros(n,s,K);
y = zeros(n,1);

fs = zeros(r,K);
cs = zeros(r,K);
Hs = zeros(s,r,K);
kappa = 20;
for kk = 1:K
    
    [fs(:,kk), cs(:,kk), Hs(:,:,kk), ~, X_gt(:,:,kk)] = get_X_with_sep_fixed_condition_number(r, s, n, kappa);
    Bs(:,:,kk) = -sqrt(3)+2*sqrt(3).*rand(n,s);
    y = y + diag(Bs(:,:,kk)*X_gt(:,:,kk));
end

%% data matrix recovery
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
maxiter = 100;
tol_rec = 1e-4;
tol_gm  = 1e-8;
tol_obj = 1e-5;


is_show = 1;

step_size = 0.6/K; %1 常数步长; 0 线搜索步长；
[~,Xs_rec, scaled_gd_recover_err] = solver_scaled_gd(y, Bs, Ls_init, Rs_init, X_gt, s, n, r, K, maxiter, tol_rec,tol_gm, tol_obj, step_size, is_show);


fprintf('Recovery Error=%f\n', scaled_gd_recover_err(end));
%% recovery frequencies via Smoothed MUSIC

% get the hankel matrix of X0;

sep = 1/n;%min(min(abs(fs_sort(2:end)-fs_sort(1:end-1))),fs_sort(1)+1-fs_sort(end));
%fprintf('sep: %.5f, RL: %.5f\n', sep, 1/n)

grid_size = sep * 0.001;
fs_recs = zeros(r,K);
for kk = 1:K
    
    X = Xs_rec(:,:,kk);
    HXn = G(X*diag(D));
    [fs_rec, ~, ~] = MUSIC(HXn.', r, grid_size);
    fs_recs(:,kk) = fs_rec';
end


%% plot frequies


figure

hold on;

fs_gt = fs(:);
fs_est=fs_recs(:);

h1 = plot(cos(2*pi*fs_gt), sin(2*pi*fs_gt),'ro', 'MarkerSize', 10);
h2 = plot(cos(2*pi*fs_est),sin(2*pi*fs_est),'rx', 'MarkerSize', 10);


t = linspace(0,1,1e4);
plot(cos(2*pi*t), sin(2*pi*t),'--b')
grid on;
box on
legend([h1, h2], {'Ground truth locations', 'Estimated locations'}, 'Location','best');

%{
myfig = gcf;
myfig.PaperUnits = 'inches';
myfig.PaperSize = [6 5.5];
myfig.PaperPosition = [0 0 6 5];
myfig.PaperPositionMode = 'manual';
figname = 'figFreqKappa20';
print( myfig, figname, '-depsc' );

%}



