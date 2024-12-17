clear all;
 
rng(1);
rand_p = rand(3,100);
rand_p(1,:) = 1+rand_p(1,:);
rand_p(2,:) = 1+2*rand_p(2,:);
rand_p(3,:) = 1+2*rand_p(3,:);
err_x = zeros(1,100);
err_y = zeros(1,100);
for i = 1:100
    [err_x(i),err_y(i)] = linear_stanford_flux(rand_p(1,i),rand_p(2,i),rand_p(3,i));
    i
end


% figure
% width = 5;     % Width in inches
% height = 4;    % Height in inches
% set(gcf,'InvertHardcopy','on');
% set(gcf,'PaperUnits', 'inches');
% papersize = get(gcf, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(gcf,'PaperPosition', myfiguresize);
% 
% hold on;
% plot(s_h_test,err_x,'o-','color',[0 0.5 0],'LineWidth',1);
% plot(s_h_test,err_y,'bo-','LineWidth',1);
% set(gca, 'YScale', 'log')
% xlabel('$p$','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% ylabel('$L_2$ error','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% legend({'$Q_1$','$Q_2$'},...
%     'interpreter','latex','FontSize',10,'Location','best');
% legend boxoff;
% saveas(gcf,'p_vs_err_flux','epsc'); 
