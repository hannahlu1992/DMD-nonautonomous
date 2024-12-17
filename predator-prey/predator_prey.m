clear all;
close all;


dt = 0.1;
%% training data
p1 = [0,2.5,5];p2 = p1;p3 = p1;
N_par = 3;
p = zeros(length(p1),length(p2),length(p3),N_par);
[p(:,:,:,1),p(:,:,:,2),p(:,:,:,3)] = ndgrid(p1,p2,p3);
size_p = size(p);
size_p1 = size(p(:,:,:,1));

%% test data reference
t = 0:dt:100;
N_T = length(t);
x0 = [3;2];
[t,w_test_true] = ode45(@(t,x) my_ode(t,x),t,x0);
w_test_true = w_test_true';

rng(1);
x0_100 = 5.*rand([2,100]);
w_test_true100 = zeros(2,N_T,100);
for m = 1:100
    [t,w_test_soln] = ode45(@(t,x) my_ode(t,x),t,x0_100(:,m));
    w_test_true100(:,:,m) = w_test_soln';
end



%% test data from modified system
p_test = zeros(N_T-1,3);
for n = 1:N_T-1
    p_test(n,1) = sin(t(n)/3)+cos(t(n))+2;
    p_test(n,2) = sin((t(n)+dt/2)/3)+cos(t(n)+dt/2)+2;
    p_test(n,3) = sin(t(n+1)/3)+cos(t(n+1))+2;
end
w_test_modified =0*w_test_true;
w_test_modified(:,1) = x0;
for n = 1:N_T-1
    [~,y] = ode45(@(t,x) my_ode_modified(t,x,p_test(n,1),p_test(n,2),...
        p_test(n,3),dt),[0,dt],w_test_modified(:,n));
    w_test_modified(:,n+1) = y(end,:)';
end

% w_test_modified100  = w_test_true100;
% for m = 1:100
%     for n = 1:N_T-1
%         [~,y] = ode45(@(t,x) my_ode_modified(t,x,p_test(n,1),p_test(n,2),...
%         p_test(n,3),dt),[0,dt],w_test_modified100(:,n,m));
%         w_test_modified100(:,n+1,m) = y(end,:)';
%     end
% end




%% training data
[N_w,~] = size(w_test_modified);
N_q = 9;
N_snap = 10;
w_train_in = zeros([N_w,N_snap,size_p1]);
w_train_out = w_train_in;
ROB_train_dmd = zeros([9,N_q,size_p1]);
Kr_train_dmd = zeros([N_q,N_q,size_p1]);
Br_train_dmd = zeros([N_q,size_p1]);

for i1 = 1:size_p1(1)
    for i2 = 1:size_p1(2)
        for i3 = 1:size_p1(3)
            rng(1);
            w_train_in(:,:,i1,i2,i3) = 5.*rand([N_w,N_snap]);
            for n = 1:N_snap
                [~,y] = ode45(@(t,x) my_ode_modified(t,x,p1(i1),p2(i2),...
                p3(i3),dt),[0,dt],w_train_in(:,n,i1,i2,i3));
                w_train_out(:,n,i1,i2,i3)= y(end,:)';
            end
        end
    end
end

for i1 = 1:size_p1(1)
    for i2 = 1:size_p1(2)
        for i3 = 1:size_p1(3)
            [ROB_train_dmd(:,:,i1,i2,i3),...
                Kr_train_dmd(:,:,i1,i2,i3),...
                Br_train_dmd(:,i1,i2,i3)] = DMD(observ(w_train_in(:,:,i1,i2,i3)),...
                observ(w_train_out(:,:,i1,i2,i3)),N_q);
%             y = ROB_train_dmd(:,:,i1,i2,i3)*(Kr_train_dmd(:,:,i1,i2,i3)...
%                 *(ROB_train_dmd(:,:,i1,i2,i3)'*observ(w_train_in(:,:,i1,i2,i3)))...
%                 +Br_train_dmd(:,i1,i2,i3));
%             y(1:2,:)-w_train_out(:,:,i1,i2,i3)
        end
    end
end
        

%% interpolate dmd PROMs
i0 = [1,1,1]; % reference point
Gamma_train_dmd = compute_Gamma_train(ROB_train_dmd,i0,p);
wr_test_dmd =0*w_test_true;wr_test_dmd(:,1) = w_test_true(:,1);

for n = 1:N_T-1
    ROB_test_interp_dmd = interpolate_ROB(Gamma_train_dmd,ROB_train_dmd,i0,p,...
        p_test(n,1),p_test(n,2),p_test(n,3),'linear');

    [ROB_test_interp_general_dmd,Kr_test_interp_dmd,Br_test_interp_dmd] = interpolate_PROM...
    (ROB_train_dmd,ROB_test_interp_dmd,Kr_train_dmd,Br_train_dmd,i0,p,...
    p_test(n,1),p_test(n,2),p_test(n,3),'spline');

    y = ROB_test_interp_general_dmd*...
        (Kr_test_interp_dmd*(ROB_test_interp_general_dmd'*observ(wr_test_dmd(:,n)))...
        +Br_test_interp_dmd);

    wr_test_dmd(:,n+1) = y(1:2);
end

wr_test_dmd100 =0*w_test_true100;wr_test_dmd100(:,1,:) = w_test_true100(:,1,:);
for m = 1:100
    for n = 1:N_T-1
        y = ROB_test_interp_general_dmd*...
            (Kr_test_interp_dmd*(ROB_test_interp_general_dmd'*observ(wr_test_dmd100(:,n,m)))...
            +Br_test_interp_dmd);
        wr_test_dmd100(:,n+1,m) = y(1:2);
    end
end
error_mean = mean(abs((w_test_true100-wr_test_dmd100)./w_test_true100),3);

% hold on;
% plot(t,error_mean(1,:));
% set(gca, 'YScale', 'log');


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
plot(t,w_test_true(1,:),'r','LineWidth',1);
plot(t,wr_test_dmd(1,:),'b--','LineWidth',1);
xlabel('$t$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('$S_1$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
legend({'reference','DRIPS'},'interpreter','latex','FontSize',10,'Location','best');
legend boxoff;

subplot(1,2,2)
hold on;
plot(t,vecnorm(w_test_true(1,:)-wr_test_dmd(1,:),2,1)./vecnorm(w_test_true(1,:),2,1),'b','LineWidth',1);
set(gca, 'YScale', 'log');
xlabel('$t$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('error','FontUnits','points','interpreter','latex',...
     'FontSize',10);
%legend({'$\Delta t = 0.05$','$\Delta t = 0.1$','$\Delta t = 0.15$'},'interpreter','latex','FontSize',10,'Location','southwest');
%legend boxoff;

saveas(gcf,'predator-prey1','epsc'); 


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
% plot(t,w_test_true(2,:),'r','LineWidth',1);
% plot(t,wr_test_dmd(2,:),'b--','LineWidth',1);
% xlabel('$t$','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% ylabel('$S_2$','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% legend({'reference','DMD'},'interpreter','latex','FontSize',10,'Location','best');
% legend boxoff;
% 
% subplot(1,2,2)
% hold on;
% plot(t,abs(w_test_true(1,:)-wr_test_dmd(1,:)),'b','LineWidth',1);
% set(gca, 'YScale', 'log');
% xlabel('$t$','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% ylabel('error','FontUnits','points','interpreter','latex',...
%      'FontSize',10);
% %legend({'$\Delta t = 0.05$','$\Delta t = 0.1$','$\Delta t = 0.15$'},'interpreter','latex','FontSize',10,'Location','southwest');
% % legend boxoff;
% 
% saveas(gcf,'predator-prey2','epsc'); 











function dx = my_ode(t,x)
dx = [x(1)-x(1)*x(2)+sin(t/3)+cos(t)+2;-x(2)+x(1)*x(2)];
end


function dx = my_ode_modified(t,x,p1,p2,p3,dt)
l1 = (t-0.5*dt)*(t-dt)/((0-0.5*dt)*(0-dt));
l2 = (t-0*dt)*(t-dt)/((0.5*dt-0*dt)*(0.5*dt-dt));
l3 = (t-0*dt)*(t-0.5*dt)/((dt-0*dt)*(dt-0.5*dt));

mu_t = p1*l1+p2*l2+p3*l3;
dx = [x(1)-x(1)*x(2)+mu_t;-x(2)+x(1)*x(2)];
end

function y = observ(x)
y = [x(1,:);x(2,:);x(1,:).^2;x(1,:).*x(2,:);x(2,:).^2;...
    x(1,:).^3;x(1,:).^2.*x(2,:);x(1,:).*x(2,:).^2;x(2,:).^3];
end


function [ROB, Kr,Br] = DMD(X1,X2,r)
    X1_tilde = [X1;ones(1,size(X1,2))];
    [V,Sigma,Z] = svd(X1_tilde,'econ');
    V = V(:,1:r+1); Sigma = Sigma(1:r+1,1:r+1); Z = Z(:,1:r+1);
    A_tilde = X2*Z/Sigma*V';
    K = A_tilde(:,1:end-1);
    B = A_tilde(:,end);
    [ROB,~,~] = svd(X1,'econ');
    ROB = ROB(:,1:r);
    Kr = ROB'*K*ROB;
    Br = ROB'*B;
end

function Gamma_train = compute_Gamma_train(ROB_train,i0,p)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,1));
    Gamma_train = zeros([N_w,N_q,size_p1]);
    matrix_pre1 = (eye(N_w)-ROB_train(:,:,i0(1),i0(2),i0(3))*ROB_train(:,:,i0(1),i0(2),i0(3))');
    for i1 = 1:size_p1(1)
        for i2 = 1:size_p1(2)
            for i3 = 1:size_p1(3)
                [U,Sigma,Z] = svd(matrix_pre1*...
            ROB_train(:,:,i1,i2,i3)*inv(ROB_train(:,:,i0(1),i0(2),i0(3))'*ROB_train(:,:,i1,i2,i3)),'econ');
                Gamma_train(:,:,i1,i2,i3) = U*diag(atan(diag(Sigma)))*Z';
            end
        end
    end
end

function ROB_test = interpolate_ROB(Gamma_train,ROB_train,i0,p,p1,p2,p3,interpolate_method)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,1));
    Gamma_test = zeros(N_w,N_q);
    parfor k = 1:N_w
        for s = 1:N_q
            Gamma_test(k,s) = interpn(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3),...
                reshape(Gamma_train(k,s,:,:,:),size_p1),p1,p2,p3,interpolate_method);
        end
    end
    [U,Sigma,Z] = svd(Gamma_test,'econ');
    ROB_test = ROB_train(:,:,i0(1),i0(2),i0(3))*Z*diag(cos(diag(Sigma)))+U*diag(sin(diag(Sigma)));
end

function [ROB_test,Kr_test,Br_test] = interpolate_PROM(ROB_train,ROB_test,Kr_train,Br_train,i0,p,p1,p2,p3,interpolate_method)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,1));
    %% Step A
    for i1 = 1:size_p1(1)
        for i2 = 1:size_p1(2)
            for i3 = 1:size_p1(3)
                P = ROB_train(:,:,i1,i2,i3)'* ROB_train(:,:,i0(1),i0(2),i0(3));
                [U,~,Z] = svd(P,'econ');
                Q = U*Z';
                Kr_train(:,:,i1,i2,i3) = Q'*Kr_train(:,:,i1,i2,i3)*Q;
                Br_train(:,i1,i2,i3) = Q'*Br_train(:,i1,i2,i3);
            end
        end
    end
    P = ROB_test'* ROB_train(:,:,i0(1),i0(2),i0(3));
    [U,~,Z] = svd(P,'econ');
    Q = U*Z';
    ROB_test = ROB_test*Q;

    %% Step B
    Gamma_K_train = manifold_log(Kr_train(:,:,i0(1),i0(2),i0(3)),Kr_train,'real');
    Gamma_K_test = zeros(N_q,N_q);
    for k = 1:N_q
        for s = 1:N_q
            Gamma_K_test(k,s) = interpn(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3),...
                reshape(Gamma_K_train(k,s,:,:,:),size_p1),p1,p2,p3,interpolate_method);
        end
    end
    Kr_test = manifold_exp(Kr_train(:,:,i0(1),i0(2),i0(3)),Gamma_K_test,'real');

    Gamma_B_train = manifold_log(Br_train(:,i0(1),i0(2),i0(3)),Br_train,'real');
    Gamma_bias_test = zeros(N_q,1);
    for k = 1:N_q
        Gamma_bias_test(k) = interpn(p(:,:,:,1),p(:,:,:,2),p(:,:,:,3),...
                reshape(Gamma_B_train(k,:,:,:),size_p1),p1,p2,p3,interpolate_method);
    end
    Br_test = manifold_exp(Br_train(:,i0(1),i0(2),i0(3)), Gamma_bias_test, 'real');
end

function Gamma = manifold_log(X,Y,matrix_type)
if strcmp(matrix_type,'real')
    Gamma = Y-X;
end
if strcmp(matrix_type,'nonsingular')
    Gamma = log(Y*inv(X));
end
if strcmp(matrix_type,'SPD')
    R = chol(X);
    Gamma = log(inv(R)*Y*inv(R));
end
end

function Y = manifold_exp(X,Gamma,matrix_type)
if strcmp(matrix_type,'real')
    Y = X+Gamma;
end
if strcmp(matrix_type,'nonsingular')
    Y = exp(Gamma)*X;
end
if strcmp(matrix_type,'SPD')
    R = chol(X);
    Y = R*exp(Gamma)*R;
end
end






