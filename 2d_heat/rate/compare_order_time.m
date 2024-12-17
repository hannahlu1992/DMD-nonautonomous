clear all;
fprintf('m=10:\n');
[e_train_rate_dmd10,e_rate_dmd10,t_offline10,t_online10] = linear_stanford_rate(10);
fprintf('m=9:\n');
[e_train_rate_dmd9,e_rate_dmd9,t_offline9,t_online9] = linear_stanford_rate(9);
fprintf('m=8:\n');
[e_train_rate_dmd8,e_rate_dmd8,t_offline8,t_online8] = linear_stanford_rate(8);
fprintf('m=7:\n');
[e_train_rate_dmd7,e_rate_dmd7,t_offline7,t_online7] = linear_stanford_rate(7);
fprintf('m=6:\n');
[e_train_rate_dmd6,e_rate_dmd6,t_offline6,t_online6] = linear_stanford_rate(6);
fprintf('m=5:\n');
[e_train_rate_dmd5,e_rate_dmd5,t_offline5,t_online5] = linear_stanford_rate(5);
fprintf('m=4:\n');
[e_train_rate_dmd4,e_rate_dmd4,t_offline4,t_online4] = linear_stanford_rate(4);
fprintf('m=3:\n');
[e_train_rate_dmd3,e_rate_dmd3,t_offline3,t_online3] = linear_stanford_rate(3);
fprintf('m=2:\n');
[e_train_rate_dmd2,e_rate_dmd2,t_offline2,t_online2] = linear_stanford_rate(2);
fprintf('m=1:\n');
[e_train_rate_dmd1,e_rate_dmd1,t_offline1,t_online1] = linear_stanford_rate(1);


t_data = 333.514541;

e_test = [e_rate_dmd1,e_rate_dmd2,e_rate_dmd3,e_rate_dmd4,e_rate_dmd5,...
    e_rate_dmd6,e_rate_dmd7,e_rate_dmd8,e_rate_dmd9,e_rate_dmd10];
e_train = [mean(reshape(e_train_rate_dmd1,[27,1])),mean(reshape(e_train_rate_dmd2,[27,1])),mean(reshape(e_train_rate_dmd3,[27,1])),...
    mean(reshape(e_train_rate_dmd4,[27,1])),mean(reshape(e_train_rate_dmd5,[27,1])),mean(reshape(e_train_rate_dmd6,[27,1])),...
    mean(reshape(e_train_rate_dmd7,[27,1])),mean(reshape(e_train_rate_dmd8,[27,1])),mean(reshape(e_train_rate_dmd9,[27,1])),...
    mean(reshape(e_train_rate_dmd10,[27,1]))];

figure
width = 5;     % Width in inches
height = 4;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

hold on;
plot(1:10,log(e_train),'bd-','LineWidth',1);
plot(1:10,log(e_test),'ro-','LineWidth',1);
xlim([1,10]);
xlabel('$m$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('average $\log(\mathcal E_{DMD})$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
legend({'train','test'},'FontUnits','points','interpreter',...
    'latex','FontSize',10,'Location','northeast');
legend boxoff;
saveas(gcf,'order_compare','epsc'); 

t_offline = [t_offline1,t_offline2,t_offline3,t_offline4,t_offline5,t_offline6,...
    t_offline7,t_offline8,t_offline9,t_offline10];
t_offline = t_offline +t_data;
t_online = [t_online1,t_online2,t_online3,t_online4,t_online5,t_online6,...
    t_online7,t_online8,t_online9,t_online10];
t_total = t_offline+t_online;


figure
width = 5;     % Width in inches
height = 4;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);


hold on;
yyaxis left
plot(1:10,t_data*ones(1,10),'b-.','LineWidth',1);
plot(1:10,t_offline,'bd-','LineWidth',1);
xlim([1,10]);

yyaxis right
plot(1:10,t_online,'ro-','LineWidth',1);
xlim([1,10]);

yyaxis left
xlabel('$m$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('computational time (s)','FontUnits','points','interpreter','latex',...
     'FontSize',10);
 
yyaxis right
ylabel('computational time (s)','FontUnits','points','interpreter','latex',...
     'FontSize',10);

ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';

legend({'generating training data','offline','online'},'FontUnits','points','interpreter',...
    'latex','FontSize',10,'Location','northwest');
legend boxoff;
saveas(gcf,'time_compare','epsc'); 
















