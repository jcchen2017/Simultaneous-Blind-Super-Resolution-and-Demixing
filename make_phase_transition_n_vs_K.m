% demixing and supre-solution
clear all;close all;clc
load('empir_probs_scaledGD_n80_160_vs_K_1_8.mat')

% with separation
colormap(gray);
imagesc(empir_probs_scaledGD_n80_160_vs_K_1_8);
fontsz = 16;
set(gca,'YDir','normal','Fontsize',fontsz) % y轴刻度设置为从下往上一次增大
%set(gca, 'ytick',80:10:240, ...
%    'yticklabel', {80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 },...
%    'xtick',[2 3 4 5 6 7 8], ...
%    'xticklabel', {2 3 4 5 6 7 8},'xticklabelmode','auto','yticklabelmode','auto');
%axis([2 8 1 17])
xticklabels({2 3 4 5 6 7 8});
%yticklabels({80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240});
%ylim([80 240]); % 设置纵坐标范围
yticks(1:4:17); % 设置纵坐标刻度

% 可选：设置纵坐标刻度标签（如果需要自定义显示）
yticklabels([80 120 160 200 240]);

cbar = colorbar;
%set(cbar,'Ticks',[0:0.2:1])
title('Scaled-GD with Separation','interpreter','latex','fontsize', fontsz)
ylabel('Number of measurement: $n$','interpreter','latex','fontsize', fontsz)
xlabel('Number of signals: $K$','interpreter','latex','fontsize', fontsz)
hold on
x = 0:0.5:8;
y =3*x-2;
plot(x,y,'red','Linewidth', 2);

myfig = gcf;
myfig.PaperUnits = 'inches';
myfig.PaperSize = [6 5.5];
myfig.PaperPosition = [0 0 6 5];
myfig.PaperPositionMode = 'manual';
figname = 'phaseTransitionScaledGDwithSep_n_vs_K';
print( myfig, figname, '-depsc' );
%}

