% clear all;
% close all;
% dt1 = 0.05;
% dt2 = 0.1;
% dt3 = 0.15;
% 
% t1 = 0:dt1:100;
% t2 = 0:dt2:100;
% t3 = 0:dt3:100;
% 
% [w_test_true1,w_test_modified1,wr_test_dmd1] = scalar(dt1);
% [w_test_true2,w_test_modified2,wr_test_dmd2] = scalar(dt2);
% [w_test_true3,w_test_modified3,wr_test_dmd3] = scalar(dt3);
% 
% 
% % w_test_true1 = interp1(t,w_test_true, t1);
% % w_test_true2 = interp1(t,w_test_true, t2);
% % w_test_true3 = interp1(t,w_test_true, t3);
% 


figure
width = 8;     % Width in inches
height = 3;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

subplot(1,2,1)
hold on;
plot(t2,w_test_true2,'r','LineWidth',1);
plot(t2,wr_test_dmd2,'b--','LineWidth',1);
xlabel('$t$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('$S$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
legend({'reference','DRIPS'},'interpreter','latex','FontSize',10,'Location','southwest');
legend boxoff;

subplot(1,2,2)
hold on;
plot(t1,vecnorm(w_test_true1-wr_test_dmd1,2,1)./(vecnorm(w_test_true1,2,1)+1e-2),'b','LineWidth',1);
plot(t2,vecnorm(w_test_true2-wr_test_dmd2,2,1)./(vecnorm(w_test_true2,2,1)+1e-2),'r--','LineWidth',1);
plot(t3,vecnorm(w_test_true3-wr_test_dmd3,2,1)./(vecnorm(w_test_true3,2,1)+1e-2),'k-.','LineWidth',1);
set(gca, 'YScale', 'log');
ylim([1e-7,1e0]);
xlabel('$t$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('error','FontUnits','points','interpreter','latex',...
     'FontSize',10);
legend({'$\Delta t = 0.05$','$\Delta t = 0.1$','$\Delta t = 0.15$'},'interpreter','latex','FontSize',10,'Location','southeast');
legend boxoff;

saveas(gcf,'scalar','epsc'); 

% figure
% width = 8;     % Width in inches
% height = 3;    % Height in inches
% set(gcf,'InvertHardcopy','on');
% set(gcf,'PaperUnits', 'inches');
% papersize = get(gcf, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(gcf,'PaperPosition', myfiguresize);
% 
% subplot(1,2,1)
% hold on;
% plot(t1,abs(w_test_true1-w_test_modified1),'b','LineWidth',1);
% plot(t2,abs(w_test_true2-w_test_modified2),'r--','LineWidth',1);
% plot(t3,abs(w_test_true3-w_test_modified3),'k-.','LineWidth',1);
% set(gca, 'YScale', 'log');
% ylim([1e-9,1e-1]);
% xlabel('$t$','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% ylabel('error','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% legend({'$\Delta t = 0.05$','$\Delta t = 0.1$','$\Delta t = 0.15$'},'interpreter','latex','FontSize',10,'Location','southwest');
% legend boxoff;
% 
% subplot(1,2,2)
% hold on;
% plot(t1,abs(w_test_modified1-wr_test_dmd1),'b','LineWidth',1);
% plot(t2,abs(w_test_modified2-wr_test_dmd2),'r--','LineWidth',1);
% plot(t3,abs(w_test_modified3-wr_test_dmd3),'k-.','LineWidth',1);
% set(gca, 'YScale', 'log');
% ylim([1e-7,1e-1]);
% xlabel('$t$','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% ylabel('error','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% legend({'$\Delta t = 0.05$','$\Delta t = 0.1$','$\Delta t = 0.15$'},'interpreter','latex','FontSize',10,'Location','southwest');
% legend boxoff;
% 
% saveas(gcf,'scalar_error','epsc'); 
