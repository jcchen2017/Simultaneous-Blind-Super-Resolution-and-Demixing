% demixing and supre-solution
clear all;close all;clc
load('data_empir_probs_scaledGD_r_vs_s_48_2024-05-17 10:37:59.mat')
load('data_empir_probs_scaledGD_r_vs_s_48_2024-05-17 10:59:19_wihtout_sep.mat')

% with separation
colormap(gray);
imagesc(empir_probs_scaledGD_r_vs_s);
fontsz = 16;
set(gca,'YDir','normal','Fontsize',fontsz) % y轴刻度设置为从下往上一次增大
%set(gca, 'ytick',[2:2:10], 'yticklabel', [2:2:10],'xtick',[2:2:18], 'xticklabel', [2:2:18],'xticklabelmode','auto','yticklabelmode','auto');
cbar = colorbar;
%set(cbar,'Ticks',[0:0.2:1])
title('Scaled-GD with Separation','interpreter','latex','fontsize', fontsz)
ylabel('Number of spikes: $r$','interpreter','latex','fontsize', fontsz)
xlabel('Dimension of subspace: $s$','interpreter','latex','fontsize', fontsz)
%hold on
%x = 0.01:0.01:20;
%y = 6./x;
%plot(x,y,'red','Linewidth', 2);
myfig = gcf;
myfig.PaperUnits = 'inches';
myfig.PaperSize = [6 5.5];
myfig.PaperPosition = [0 0 6 5];
myfig.PaperPositionMode = 'manual';
figname = 'phaseTransitionScaledGDwithSep';
print( myfig, figname, '-depsc' );



% without separation
figure()
colormap(gray);
imagesc(empir_probs_scaledGD_r_vs_s_without_sep);
fontsz = 16;
set(gca,'YDir','normal','Fontsize',fontsz) % y轴刻度设置为从下往上一次增大
%set(gca, 'ytick',[2:2:10], 'yticklabel', [2:2:10],'xtick',[2:2:18], 'xticklabel', [2:2:18],'xticklabelmode','auto','yticklabelmode','auto');
cbar = colorbar;
%set(cbar,'Ticks',[0:0.2:1])
title('Scaled-GD without Separation','interpreter','latex','fontsize', fontsz)
ylabel('Number of spikes: $r$','interpreter','latex','fontsize', fontsz)
xlabel('Dimension of subspace: $s$','interpreter','latex','fontsize', fontsz)
%hold on
%x = 0.01:0.01:20;
%y = 6./x;
%plot(x,y,'red','Linewidth', 2);
myfig = gcf;
myfig.PaperUnits = 'inches';
myfig.PaperSize = [6 5.5];
myfig.PaperPosition = [0 0 6 5];
myfig.PaperPositionMode = 'manual';
figname = 'phaseTransitionScaledGDwithoutSep';
print( myfig, figname, '-depsc' );

