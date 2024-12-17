clear all;

order = 7;
 
rng(1);
rand_p = rand(3,100);
rand_p(1,:) = 1+rand_p(1,:);
rand_p(2,:) = 1+2*rand_p(2,:);
rand_p(3,:) = 1+2*rand_p(3,:);
e_dmd = zeros(1,100);
for i = 1:100
    e_dmd(i) = linear_stanford_rate(rand_p(1,i),rand_p(2,i),rand_p(3,i),order);
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
% plot(s_h_test,err,'ko-','LineWidth',1);
% set(gca, 'YScale', 'log')
% xlabel('$p$','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% ylabel('$L_2$ error','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% saveas(gcf,'p_vs_err_rate','epsc'); 
