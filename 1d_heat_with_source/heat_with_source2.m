clear all;
close all;


dt = 0.05;
%% training data
p1 = [0,0.5,1];p2 = p1;p3 = p1;p4 = p1;p5 = p1;p6 = [0,1]; p7 = [0.05,0.5];
N_par = 7;
p = zeros(length(p1),length(p2),length(p3),length(p4),length(p5),length(p6),length(p7),N_par);
[p(:,:,:,:,:,:,:,1),p(:,:,:,:,:,:,:,2),p(:,:,:,:,:,:,:,3),p(:,:,:,:,:,:,:,4),...
    p(:,:,:,:,:,:,:,5),p(:,:,:,:,:,:,:,6),p(:,:,:,:,:,:,:,7)] = ndgrid(p1,p2,p3,p4,p5,p6,p7);
size_p = size(p);
size_p1 = size(p(:,:,:,:,:,:,:,1));

%% test data reference
xmesh = linspace(0,1,22);
t = 0:dt:8;N_T = length(t);

m = 0;
mu = 1;
sigma = 0.5;
w_test_true = pdepe(m,@(x,t,u,dudx) heatpde(x,t,u,dudx,mu,sigma),@heatic,@heatbc,xmesh,t);
w_test_true = w_test_true';

%% test data from modified system
p_test = zeros(N_T-1,N_par);
for n = 1:N_T-1
    p_test(n,1) = my_fun(t(n));
    p_test(n,2) = my_fun(t(n)+dt/4);
    p_test(n,3) = my_fun(t(n)+dt/2);
    p_test(n,4) = my_fun(t(n)+dt/4*3);
    p_test(n,5) = my_fun(t(n+1));
    p_test(n,6) = mu;
    p_test(n,7) = sigma;
end

w_test_modified = 0*w_test_true;
w_test_modified(:,1) = w_test_true(:,1);
tau = 0:dt/2:dt;
for n = 1:N_T-1
    icarg = @(x)interp1(xmesh,w_test_modified(:,n)',x);
    y = pdepe(m,@(x,tau,u,dudx) heatpde_modified(x,tau,u,dudx,mu,sigma,p_test(n,1),...
        p_test(n,2),p_test(n,3),p_test(n,4),p_test(n,5),dt),@(x) heatic_modified(x,icarg),@heatbc,xmesh,tau);
    w_test_modified(:,n+1) = y(end,:)';
end

figure
subplot(1,3,1)
imagesc(w_test_true)
colorbar

subplot(1,3,2)
imagesc(w_test_modified)
colorbar

subplot(1,3,3)
imagesc(abs(w_test_true-w_test_modified))
colorbar

%% training data
[N_w,~] = size(w_test_modified);
N_q = N_w;
N_snap = 25;
w_train_in = zeros([N_w,N_snap,size_p1]);
w_train_out = w_train_in;
ROB_train_dmd = zeros([N_w,N_q,size_p1]);
Kr_train_dmd = zeros([N_q,N_q,size_p1]);
Br_train_dmd = zeros([N_q,size_p1]);

for i1 = 1:size_p1(1)
    for i2 = 1:size_p1(2)
        for i3 = 1:size_p1(3)
            for i4 = 1:size_p1(4)
                for i5 = 1:size_p1(5)
                    for i6 = 1:size_p1(6)
                        for i7 = 1:size_p1(7)
                            w_train_in(:,:,i1,i2,i3,i4,i5,i6,i7) = 2.*rand([N_w,N_snap]);
                            for n = 1:N_snap
                                icarg = @(x)interp1(xmesh,w_train_in(:,n,i1,i2,i3,i4,i5,i6,i7)',x);
                                y = pdepe(m,@(x,tau,u,dudx) heatpde_modified(x,tau,u,dudx,p6(i6),p7(i7),p1(i1),p2(i2),...
                                p3(i3),p4(i4),p5(i5),dt),@(x) heatic_modified(x,icarg),@heatbc,xmesh,tau);
                                w_train_out(:,n,i1,i2,i3,i4,i5,i6,i7)= y(end,:)';
                            end
                        end
                    end
                end
            end
        end
    end
end

for i1 = 1:size_p1(1)
    for i2 = 1:size_p1(2)
        for i3 = 1:size_p1(3)
            for i4 = 1:size_p1(4)
                for i5 = 1:size_p1(5)
                    for i6 = 1:size_p1(6)
                        for i7 = 1:size_p1(7)
                            [ROB_train_dmd(:,:,i1,i2,i3,i4,i5,i6,i7),...
                                Kr_train_dmd(:,:,i1,i2,i3,i4,i5,i6,i7),...
                                Br_train_dmd(:,i1,i2,i3,i4,i5,i6,i7)] = DMD(w_train_in(:,:,i1,i2,i3,i4,i5,i6,i7),...
                                w_train_out(:,:,i1,i2,i3,i4,i5,i6,i7),N_q);
                        end
                    end
                end
            end
        end
    end
end
        
%% interpolate dmd PROMs
i0 = [1,1,1,1,1,1,1]; % reference point
Gamma_train_dmd = compute_Gamma_train(ROB_train_dmd,i0,p);
wr_test_dmd =0*w_test_true;wr_test_dmd(:,1) = w_test_true(:,1);

for n = 1:N_T-1
    ROB_test_interp_dmd = interpolate_ROB(Gamma_train_dmd,ROB_train_dmd,i0,p,...
        p_test(n,1),p_test(n,2),p_test(n,3),p_test(n,4),p_test(n,5),p_test(n,6),p_test(n,7),'linear');
    [ROB_test_interp_general_dmd,Kr_test_interp_dmd,Br_test_interp_dmd] = interpolate_PROM...
    (ROB_train_dmd,ROB_test_interp_dmd,Kr_train_dmd,Br_train_dmd,i0,p,...
    p_test(n,1),p_test(n,2),p_test(n,3),p_test(n,4),p_test(n,5),p_test(n,6),p_test(n,7),'spline');

    wr_test_dmd(:,n+1) = ROB_test_interp_general_dmd*...
        (Kr_test_interp_dmd*(ROB_test_interp_general_dmd'*wr_test_dmd(:,n))...
        +Br_test_interp_dmd);
end

figure
width = 8;     % Width in inches
height = 6;    % Height in inches
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

subplot(2,2,1)
hold on;
plot(t,w_test_true(11,:),'r','LineWidth',1);
plot(t,wr_test_dmd(11,:),'b.','Markersize',10);
xlabel('$t$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('$S$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
legend({'reference','DRIPS'},'interpreter','latex','FontSize',10,'Location','northwest');
legend boxoff;

subplot(2,2,2)
hold on;
plot(xmesh,w_test_true(:,11),'r','LineWidth',1);
plot(xmesh,wr_test_dmd(:,11),'b.','Markersize',10);
xlabel('$x$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('$S$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
legend({'reference','DRIPS'},'interpreter','latex','FontSize',10,'Location','northwest');
legend boxoff;

subplot(2,2,3)
[~,c] = contour(xmesh,t,w_test_true',20);
c.LineWidth = 1;
xlabel('$x$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('$t$','FontUnits','points','interpreter','latex',...
     'FontSize',10);

subplot(2,2,4)
[~,c] = contour(xmesh,t,wr_test_dmd',20);
c.LineWidth = 1;
xlabel('$x$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('$t$','FontUnits','points','interpreter','latex',...
     'FontSize',10);

saveas(gcf,'heat_with_source','epsc'); 


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

hold on;
plot(t,vecnorm(w_test_true-wr_test_dmd,2,1)./(vecnorm(w_test_true,2,1)+1e-2),'k','LineWidth',1);
xlabel('$t$','FontUnits','points','interpreter','latex',...
     'FontSize',10);
ylabel('error','FontUnits','points','interpreter','latex',...
     'FontSize',10);
set(gca, 'YScale', 'log')
saveas(gcf,'heat_with_source_err','epsc'); 

function alpha = my_fun(t)
alpha = t-floor(t);
end



function [c,f,s] = heatpde(x,t,u,dudx,mu,sigma)
c = 1;
f = dudx;
alpha = my_fun(t);
s = alpha*exp(-(x-mu)^2/sigma^2);
end

function u0 = heatic(x)
u0 = 0.15*exp(x/10);
end

function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur;
qr = 0;
end

function [c,f,s] = heatpde_modified(x,t,u,dudx,mu,sigma,p1,p2,p3,p4,p5,dt)
c = 1;
f = dudx;

l1 = (t-0.25*dt)*(t-0.5*dt)*(t-0.75*dt)*(t-dt)/((0-0.25*dt)*(0-0.5*dt)*(0-0.75*dt)*(0-dt));
l2 = (t-0*dt)*(t-0.5*dt)*(t-0.75*dt)*(t-dt)/((0.25*dt-0*dt)*(0.25*dt-0.5*dt)*(0.25*dt-0.75*dt)*(0.25*dt-dt));
l3 = (t-0*dt)*(t-0.25*dt)*(t-0.75*dt)*(t-dt)/((0.5*dt-0*dt)*(0.5*dt-0.25*dt)*(0.5*dt-0.75*dt)*(0.5*dt-dt));
l4 = (t-0*dt)*(t-0.25*dt)*(t-0.5*dt)*(t-dt)/((0.75*dt-0*dt)*(0.75*dt-0.25*dt)*(0.75*dt-0.5*dt)*(0.75*dt-dt));
l5 = (t-0*dt)*(t-0.25*dt)*(t-0.5*dt)*(t-0.75*dt)/((dt-0*dt)*(dt-0.25*dt)*(dt-0.5*dt)*(dt-0.75*dt));


alpha_t = p1*l1+p2*l2+p3*l3+p4*l4+p5*l5; 
s = alpha_t*exp(-(x-mu)^2/sigma^2);
end

function u0 = heatic_modified(x,icarg)
u0 = icarg(x);
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
    size_p1 = size(p(:,:,:,:,:,:,:,1));
    Gamma_train = zeros([N_w,N_q,size_p1]);
    matrix_pre1 = (eye(N_w)-ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5),i0(6),i0(7))...
        *ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5),i0(6),i0(7))');
    for i1 = 1:size_p1(1)
        for i2 = 1:size_p1(2)
            for i3 = 1:size_p1(3)
                for i4 = 1:size_p1(4)
                    for i5 = 1:size_p1(5)
                        for i6 = 1:size_p1(6)
                            for i7 = 1:size_p1(7)
                [U,Sigma,Z] = svd(matrix_pre1*...
            ROB_train(:,:,i1,i2,i3,i4,i5,i6,i7)*inv(ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5),i0(6),i0(7))'*ROB_train(:,:,i1,i2,i3,i4,i5,i6,i7)),'econ');
                Gamma_train(:,:,i1,i2,i3,i4,i5,i6,i7) = U*diag(atan(diag(Sigma)))*Z';
                            end
                        end
                    end
                end
            end
        end
    end
end

function ROB_test = interpolate_ROB(Gamma_train,ROB_train,i0,p,p1,p2,p3,p4,p5,p6,p7,interpolate_method)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,:,:,:,:,1));
    Gamma_test = zeros(N_w,N_q);
    parfor k = 1:N_w
        for s = 1:N_q
            Gamma_test(k,s) = interpn(p(:,:,:,:,:,:,:,1),p(:,:,:,:,:,:,:,2),p(:,:,:,:,:,:,:,3),...
                p(:,:,:,:,:,:,:,4),p(:,:,:,:,:,:,:,5),p(:,:,:,:,:,:,:,6),p(:,:,:,:,:,:,:,7),...
                reshape(Gamma_train(k,s,:,:,:,:,:,:,:),size_p1),p1,p2,p3,p4,p5,p6,p7,interpolate_method);
        end
    end
    [U,Sigma,Z] = svd(Gamma_test,'econ');
    ROB_test = ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5),i0(6),i0(7))*Z*diag(cos(diag(Sigma)))+U*diag(sin(diag(Sigma)));
end

function [ROB_test,Kr_test,Br_test] = interpolate_PROM(ROB_train,ROB_test,Kr_train,Br_train,i0,p,p1,p2,p3,p4,p5,p6,p7,interpolate_method)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,:,:,:,:,1));
    %% Step A
    for i1 = 1:size_p1(1)
        for i2 = 1:size_p1(2)
            for i3 = 1:size_p1(3)
                for i4 = 1:size_p1(4)
                    for i5 = 1:size_p1(5)
                        for i6 = 1:size_p1(6)
                            for i7 = 1:size_p1(7)
                        P = ROB_train(:,:,i1,i2,i3,i4,i5,i6,i7)'* ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5),i0(6),i0(7));
                        [U,~,Z] = svd(P,'econ');
                        Q = U*Z';
                        Kr_train(:,:,i1,i2,i3,i4,i5,i6,i7) = Q'*Kr_train(:,:,i1,i2,i3,i4,i5,i6,i7)*Q;
                        Br_train(:,i1,i2,i3,i4,i5,i6,i7) = Q'*Br_train(:,i1,i2,i3,i4,i5,i6,i7);
                            end
                        end
                    end
                end
            end
        end
    end
    P = ROB_test'* ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5),i0(6),i0(7));
    [U,~,Z] = svd(P,'econ');
    Q = U*Z';
    ROB_test = ROB_test*Q;

    %% Step B
    Gamma_K_train = manifold_log(Kr_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5),i0(6),i0(7)),Kr_train,'real');
    Gamma_K_test = zeros(N_q,N_q);
    for k = 1:N_q
        for s = 1:N_q
            Gamma_K_test(k,s) = interpn(p(:,:,:,:,:,:,:,1),p(:,:,:,:,:,:,:,2),p(:,:,:,:,:,:,:,3),...
                p(:,:,:,:,:,:,:,4),p(:,:,:,:,:,:,:,5),p(:,:,:,:,:,:,:,6),p(:,:,:,:,:,:,:,7),...
                reshape(Gamma_K_train(k,s,:,:,:,:,:,:,:),size_p1),p1,p2,p3,p4,p5,p6,p7,interpolate_method);
        end
    end
    Kr_test = manifold_exp(Kr_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5),i0(6),i0(7)),Gamma_K_test,'real');

    Gamma_B_train = manifold_log(Br_train(:,i0(1),i0(2),i0(3),i0(4),i0(5),i0(6),i0(7)),Br_train,'real');
    Gamma_bias_test = zeros(N_q,1);
    for k = 1:N_q
        Gamma_bias_test(k) = interpn(p(:,:,:,:,:,:,:,1),p(:,:,:,:,:,:,:,2),p(:,:,:,:,:,:,:,3),...
                            p(:,:,:,:,:,:,:,4),p(:,:,:,:,:,:,:,5),p(:,:,:,:,:,:,:,6),p(:,:,:,:,:,:,:,7),...
                reshape(Gamma_B_train(k,:,:,:,:,:,:,:),size_p1),p1,p2,p3,p4,p5,p6,p7,interpolate_method);
    end
    Br_test = manifold_exp(Br_train(:,i0(1),i0(2),i0(3),i0(4),i0(5),i0(6),i0(7)), Gamma_bias_test, 'real');
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






