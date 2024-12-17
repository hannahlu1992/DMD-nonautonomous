close all;
clear all;

N_T= 500;
t = linspace(0,5000,N_T);dt  = t(2)-t(1);

%% training data
p1 = [0,1,2];p2 = p1;p3 = p1;p4 = p1;p5 = p1;
N_par = 5;
p = zeros(length(p1),length(p2),length(p3),length(p4),length(p5),N_par);
[p(:,:,:,:,:,1),p(:,:,:,:,:,2),p(:,:,:,:,:,3),p(:,:,:,:,:,4),p(:,:,:,:,:,5)] = ndgrid(p1,p2,p3,p4,p5);
size_p = size(p);
size_p1 = size(p(:,:,:,:,:,1));

%% test data reference
[w_test_true,qx_test_true,qy_test_true,rate_test_true] = stanford(t);

%% test data from modified system
p_test = zeros(N_T-1,N_par);
for n = 1:N_T-1
    p_test(n,1) = 1-cos(t(n)/800*pi);
    p_test(n,2) = 1-cos((t(n)+dt/4)/800*pi);
    p_test(n,3) = 1-cos((t(n)+dt/2)/800*pi);
    p_test(n,4) = 1-cos((t(n)+dt/4*3)/800*pi);
    p_test(n,5) = 1-cos(t(n+1)/800*pi);
end

w_test_modified = 0*w_test_true;
w_test_modified(:,1) = w_test_true(:,1);
qx_test_modified = 0*qx_test_true;
qx_test_modified(:,1) = qx_test_true(:,1);
qy_test_modified = 0*qy_test_true;
qy_test_modified(:,1) = qy_test_true(:,1);
rate_test_modified = 0*rate_test_true;
rate_test_modified(:,1) = rate_test_true(:,1);


% tau = 0:dt/2:dt;
% [w_test_modified(:,2),qx_test_modified(:,2),qy_test_modified(:,2),rate_test_modified(:,2)] = stanford_modified_step1(tau,...
%     p_test(n,1),p_test(n,2),p_test(n,3),p_test(n,4),p_test(n,5),dt);
% 
% for n = 2:N_T-1
%     n
%     [w_test_modified(:,n+1),qx_test_modified(:,n+1),qy_test_modified(:,n+1),rate_test_modified(:,n+1),~,~,~]  = stanford_modified(tau,w_test_modified(:,n),...
%     p_test(n,1),p_test(n,2),p_test(n,3),p_test(n,4),p_test(n,5),dt);
% end

load('test_data.mat');
% save('test_data.mat', 'w_test_true','w_test_modified','qx_test_true','qx_test_modified','qy_test_true','qy_test_modified','rate_test_true','rate_test_modified');
% visualize(w_test_true);
% visualize(w_test_modified);
% visualize(abs(w_test_modified-w_test_true));




%% training data
[N_w,~] = size(w_test_modified);
N_snap = N_T;
N_q = 20;
w_train = zeros([N_w,N_snap,size_p1]);
qx_train = zeros([N_w,N_snap,size_p1]);
qy_train = zeros([N_w,N_snap,size_p1]);
rate_train= zeros([1,N_snap,size_p1]);


ROB_train_dmd = zeros([N_w,N_q,size_p1]);
Kr_train_dmd = zeros([N_q,N_q,size_p1]);
Br_train_dmd = zeros([N_q,size_p1]);

% for i1 = 1:size_p1(1)
%     for i2 = 1:size_p1(2)
%         for i3 = 1:size_p1(3)
%             for i4 = 1:size_p1(4)
%                 for i5 = 1:size_p1(5)
%                     [i1,i2,i3,i4,i5]
%                     w_train(:,1,i1,i2,i3,i4,i5) = w_test_true(:,1);
%                     qx_train(:,1,i1,i2,i3,i4,i5) = qx_test_true(:,1);
%                     qy_train(:,1,i1,i2,i3,i4,i5) = qy_test_true(:,1);
%                     rate_train(:,1,i1,i2,i3,i4,i5) = rate_test_true(:,1);
%                     [w_train(:,2,i1,i2,i3,i4,i5),qx_train(:,2,i1,i2,i3,i4,i5),...
%                         qy_train(:,2,i1,i2,i3,i4,i5),rate_train(:,2,i1,i2,i3,i4,i5)] ...
%                         = stanford_modified_step1(tau,...
%                             p1(i1),p2(i2),p3(i3),p4(i4),p5(i5),dt);
%                     for n = 2:N_snap-1
%                         [w_train(:,n+1,i1,i2,i3,i4,i5),qx_train(:,n+1,i1,i2,i3,i4,i5),...
%                         qy_train(:,n+1,i1,i2,i3,i4,i5),rate_train(:,n+1,i1,i2,i3,i4,i5),~,~,~] ...
%                         = stanford_modified(tau,w_train(:,n,i1,i2,i3,i4,i5),...
%                             p1(i1),p2(i2),p3(i3),p4(i4),p5(i5),dt);
%                     end
%                 end
%             end
%         end
%     end
% end

% save('train_data.mat','w_train','qx_train','qy_train','rate_train');
load('train_data.mat');


w_train_dmd = w_train;
e_train_dmd = zeros(size_p1);

for i1 = 1:size_p1(1)
    for i2 = 1:size_p1(2)
        for i3 = 1:size_p1(3)
            for i4 = 1:size_p1(4)
                for i5 = 1:size_p1(5)
                    [ROB_train_dmd(:,:,i1,i2,i3,i4,i5), Kr_train_dmd(:,:,i1,i2,i3,i4,i5),Br_train_dmd(:,i1,i2,i3,i4,i5)] = DMD(w_train(:,1:end-1,i1,i2,i3,i4,i5),w_train(:,2:end,i1,i2,i3,i4,i5),N_q);
                    for n = 1:N_snap-1
                        w_train_dmd(:,n+1,i1,i2,i3,i4,i5) = ROB_train_dmd(:,:,i1,i2,i3,i4,i5)*...
                        (Kr_train_dmd(:,:,i1,i2,i3,i4,i5)*(ROB_train_dmd(:,:,i1,i2,i3,i4,i5)'*w_train_dmd(:,n,i1,i2,i3,i4,i5))...
                        +Br_train_dmd(:,i1,i2,i3,i4,i5));
                    end
                    e_train_dmd(i1,i2,i3,i4,i5) = norm(w_train(:,:,i1,i2,i3,i4,i5)-w_train_dmd(:,:,i1,i2,i3,i4,i5))./norm(w_train(:,:,i1,i2,i3,i4,i5));
                end
            end
        end
    end
end

%% interpolate dmd PROMs
i0 = [1,1,1,1,1]; % reference point
Gamma_train_dmd = compute_Gamma_train(ROB_train_dmd,i0,p);
wr_test_dmd = 0*w_test_true;wr_test_dmd(:,1) = w_test_true(:,1);

for n = 1:N_T-1
    n
        ROB_test_interp_dmd = interpolate_ROB(Gamma_train_dmd,ROB_train_dmd,i0,p,...
        p_test(n,1),p_test(n,2),p_test(n,3),p_test(n,4),p_test(n,5),'linear');
    [ROB_test_interp_general_dmd,Kr_test_interp_dmd,Br_test_interp_dmd] = interpolate_PROM...
    (ROB_train_dmd,ROB_test_interp_dmd,Kr_train_dmd,Br_train_dmd,i0,p,...
    p_test(n,1),p_test(n,2),p_test(n,3),p_test(n,4),p_test(n,5),'spline');

    wr_test_dmd(:,n+1) = ROB_test_interp_general_dmd*...
        (Kr_test_interp_dmd*(ROB_test_interp_general_dmd'*wr_test_dmd(:,n))...
        +Br_test_interp_dmd);
end

visualize(w_test_true);
visualize(real(wr_test_dmd));
visualize(log(abs(w_test_true-real(wr_test_dmd))));
e_dmd = norm(w_test_true-real(wr_test_dmd))./norm(w_test_true);

visualize(log(abs(w_train(:,:,1,1,1,1,1)-real(w_train_dmd(:,:,1,1,1,1,1)))));

function [T,qx,qy,w] = stanford(tlist)
    x = [167,167,349,349,462,462,260,167,167,258,561,652,652,455,455,341,341,562,652,652,547,269];
    y = [146,295,295,226,226,327,327,417,691,782,782,690,542,542,603,603,508,508,418,146,42,42];
    y = y-10;

    thermalmodel = createpde('thermal','transient');
    R1 =  [3;4;0;800;800;0;0;0;800;800];
    n = length(x);
    P1 = [2;n;x';y'];
    R1 = [R1;zeros(length(P1)-length(R1),1)];
    gd = [R1,P1]; 
    sf = 'R1-P1';
    ns = char('R1','P1');
    ns = ns';
    g = decsg(gd,sf,ns);
    geometryFromEdges(thermalmodel,g);
    thermalProperties(thermalmodel,'ThermalConductivity',1,...
                                   'MassDensity',1,...
                                   'SpecificHeat',1);

    internalHeatSource(thermalmodel,0);

    %% BC
    thermalBC(thermalmodel,'Edge',1,'Temperature',0);
    thermalBC(thermalmodel,'Edge',23,'Temperature',0);
    thermalBC(thermalmodel,'Edge',3:22,'Temperature',@my_transientBC);
    thermalBC(thermalmodel,'Edge',25:26,'Temperature',@my_transientBC);

    %% IC
    thermalIC(thermalmodel,1);

    msh = generateMesh(thermalmodel,'Hmax',50);

    pos1 = [0.08 0.2 0.19 0.7];
    pos2 = [0.41 0.2 0.19 0.7];
    pos3 = [0.74 0.2 0.19 0.7];
    c1_pos = [0.29 0.2 0.01 0.7];
    c2_pos = [0.62 0.2 0.01 0.7];
    c3_pos = [0.95 0.2 0.01 0.7];


%     thermalmodel.SolverOptions.ReportStatistics ='on';
    result = solve(thermalmodel,tlist);
    T = result.Temperature;
    [qx,qy] = evaluateHeatFlux(result);
    Q = zeros(length(tlist),26);
    for i = 1:26
        Q(:,i) = evaluateHeatRate(result,'Edge',i);
    end
    w = Q(:,3)';
end

function leftTemp = my_transientBC(~, state)
if(isnan(state.time))
  leftTemp = NaN;
else
  leftTemp = 1-cos(state.time/800*pi);
end
end

function [Tend,qxend,qyend,wend] = stanford_modified_step1(tlist,p1,p2,p3,p4,p5,dt)
    x = [167,167,349,349,462,462,260,167,167,258,561,652,652,455,455,341,341,562,652,652,547,269];
    y = [146,295,295,226,226,327,327,417,691,782,782,690,542,542,603,603,508,508,418,146,42,42];
    y = y-10;

    thermalmodel = createpde('thermal','transient');
    R1 =  [3;4;0;800;800;0;0;0;800;800];
    n = length(x);
    P1 = [2;n;x';y'];
    R1 = [R1;zeros(length(P1)-length(R1),1)];
    gd = [R1,P1]; 
    sf = 'R1-P1';
    ns = char('R1','P1');
    ns = ns';
    g = decsg(gd,sf,ns);
    geometryFromEdges(thermalmodel,g);
    thermalProperties(thermalmodel,'ThermalConductivity',1,...
                                   'MassDensity',1,...
                                   'SpecificHeat',1);

    internalHeatSource(thermalmodel,0);

    %% BC
    thermalBC(thermalmodel,'Edge',1,'Temperature',0);
    thermalBC(thermalmodel,'Edge',23,'Temperature',0);
    thermalBC(thermalmodel,'Edge',3:22,'Temperature',@(~,state) my_transientBC_modified(state,state,p1,p2,p3,p4,p5,dt));
    thermalBC(thermalmodel,'Edge',25:26,'Temperature',@(~,state) my_transientBC_modified(state,state,p1,p2,p3,p4,p5,dt));

    %% IC
    thermalIC(thermalmodel,1);

    msh = generateMesh(thermalmodel,'Hmax',50);

    pos1 = [0.08 0.2 0.19 0.7];
    pos2 = [0.41 0.2 0.19 0.7];
    pos3 = [0.74 0.2 0.19 0.7];
    c1_pos = [0.29 0.2 0.01 0.7];
    c2_pos = [0.62 0.2 0.01 0.7];
    c3_pos = [0.95 0.2 0.01 0.7];


%     thermalmodel.SolverOptions.ReportStatistics ='on';
    result = solve(thermalmodel,tlist);
    T = result.Temperature;
    Tend = T(:,end);
    [qx,qy] = evaluateHeatFlux(result);
    qxend = qx(:,end);
    qyend = qy(:,end);
    Q = zeros(length(tlist),26);
    for i = 1:26
        Q(:,i) = evaluateHeatRate(result,'Edge',i);
    end
    w = Q(:,3)';
    wend = w(:,end);
end

function [Tend,qxend,qyend,wend,qx1,qy1,w1]  = stanford_modified(tlist,Told,p1,p2,p3,p4,p5,dt)
    x = [167,167,349,349,462,462,260,167,167,258,561,652,652,455,455,341,341,562,652,652,547,269];
    y = [146,295,295,226,226,327,327,417,691,782,782,690,542,542,603,603,508,508,418,146,42,42];
    y = y-10;

    thermalmodel = createpde('thermal','transient');
    R1 =  [3;4;0;800;800;0;0;0;800;800];
    n = length(x);
    P1 = [2;n;x';y'];
    R1 = [R1;zeros(length(P1)-length(R1),1)];
    gd = [R1,P1]; 
    sf = 'R1-P1';
    ns = char('R1','P1');
    ns = ns';
    g = decsg(gd,sf,ns);
    geometryFromEdges(thermalmodel,g);
    thermalProperties(thermalmodel,'ThermalConductivity',1,...
                                   'MassDensity',1,...
                                   'SpecificHeat',1);

    internalHeatSource(thermalmodel,0);

    %% BC
    thermalBC(thermalmodel,'Edge',1,'Temperature',0);
    thermalBC(thermalmodel,'Edge',23,'Temperature',0);
    thermalBC(thermalmodel,'Edge',3:22,'Temperature',@(~,state) my_transientBC_modified(state,state,p1,p2,p3,p4,p5,dt));
    thermalBC(thermalmodel,'Edge',25:26,'Temperature',@(~,state) my_transientBC_modified(state,state,p1,p2,p3,p4,p5,dt));

    msh = generateMesh(thermalmodel,'Hmax',50);
    %% IC
    interpolant = scatteredInterpolant(thermalmodel.Mesh.Nodes(1,:)',...
        thermalmodel.Mesh.Nodes(2,:)',Told);
    thermalIC(thermalmodel,@(region) icFcn(region,interpolant));

%     thermalmodel.SolverOptions.ReportStatistics ='on';
    result = solve(thermalmodel,tlist);
    T = result.Temperature;
    Tend = T(:,end);
    [qx,qy] = evaluateHeatFlux(result);
    qxend = qx(:,end);
    qyend = qy(:,end);
    qx1 = qx(:,1);
    qy1 = qy(:,1);
    Q = zeros(length(tlist),26);
    for i = 1:26
        Q(:,i) = evaluateHeatRate(result,'Edge',i);
    end
    w = Q(:,3)';
    wend = w(:,end);
    w1 = w(:,1);
end

function u0 = icFcn(region,interpolant)
u0 = interpolant(region.x',region.y');
end

function leftTemp = my_transientBC_modified(~,state,p1,p2,p3,p4,p5,dt)
t = state.time;
l1 = (t-0.25*dt)*(t-0.5*dt)*(t-0.75*dt)*(t-dt)/((0-0.25*dt)*(0-0.5*dt)*(0-0.75*dt)*(0-dt));
l2 = (t-0*dt)*(t-0.5*dt)*(t-0.75*dt)*(t-dt)/((0.25*dt-0*dt)*(0.25*dt-0.5*dt)*(0.25*dt-0.75*dt)*(0.25*dt-dt));
l3 = (t-0*dt)*(t-0.25*dt)*(t-0.75*dt)*(t-dt)/((0.5*dt-0*dt)*(0.5*dt-0.25*dt)*(0.5*dt-0.75*dt)*(0.5*dt-dt));
l4 = (t-0*dt)*(t-0.25*dt)*(t-0.5*dt)*(t-dt)/((0.75*dt-0*dt)*(0.75*dt-0.25*dt)*(0.75*dt-0.5*dt)*(0.75*dt-dt));
l5 = (t-0*dt)*(t-0.25*dt)*(t-0.5*dt)*(t-0.75*dt)/((dt-0*dt)*(dt-0.25*dt)*(dt-0.5*dt)*(dt-0.75*dt));

if(isnan(state.time))
  leftTemp = NaN;
else
  leftTemp = p1*l1+p2*l2+p3*l3+p4*l4+p5*l5; 
end
end




function Gamma_train = compute_Gamma_train(ROB_train,i0,p)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,:,:,1));
    Gamma_train = zeros([N_w,N_q,size_p1]);
    matrix_pre1 = (eye(N_w)-ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5))...
        *ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5))');
    for i1 = 1:size_p1(1)
        for i2 = 1:size_p1(2)
            for i3 = 1:size_p1(3)
                for i4 = 1:size_p1(4)
                    for i5 = 1:size_p1(5)
                        [U,Sigma,Z] = svd(matrix_pre1*...
                    ROB_train(:,:,i1,i2,i3,i4,i5)*inv(ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5))'*ROB_train(:,:,i1,i2,i3,i4,i5)),'econ');
                        Gamma_train(:,:,i1,i2,i3,i4,i5) = U*diag(atan(diag(Sigma)))*Z';
                    end
                end
            end
        end
    end
end

function ROB_test = interpolate_ROB(Gamma_train,ROB_train,i0,p,p1,p2,p3,p4,p5,interpolate_method)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,:,:,1));
    Gamma_test = zeros(N_w,N_q);
    parfor k = 1:N_w
        for s = 1:N_q
            Gamma_test(k,s) = interpn(p(:,:,:,:,:,1),p(:,:,:,:,:,2),p(:,:,:,:,:,3),...
                p(:,:,:,:,:,4),p(:,:,:,:,:,5),...
                reshape(Gamma_train(k,s,:,:,:,:,:),size_p1),p1,p2,p3,p4,p5,interpolate_method);
        end
    end
    [U,Sigma,Z] = svd(Gamma_test,'econ');
    ROB_test = ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5))*Z*diag(cos(diag(Sigma)))+U*diag(sin(diag(Sigma)));
end

function [ROB_test,Kr_test,Br_test] = interpolate_PROM(ROB_train,ROB_test,Kr_train,Br_train,i0,p,p1,p2,p3,p4,p5,interpolate_method)
    [N_w,N_q,~] = size(ROB_train);
    size_p1 = size(p(:,:,:,:,:,1));
    %% Step A
    for i1 = 1:size_p1(1)
        for i2 = 1:size_p1(2)
            for i3 = 1:size_p1(3)
                for i4 = 1:size_p1(4)
                    for i5 = 1:size_p1(5)
                        P = ROB_train(:,:,i1,i2,i3,i4,i5)'* ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5));
                        [U,~,Z] = svd(P,'econ');
                        Q = U*Z';
                        Kr_train(:,:,i1,i2,i3,i4,i5) = Q'*Kr_train(:,:,i1,i2,i3,i4,i5)*Q;
                        Br_train(:,i1,i2,i3,i4,i5) = Q'*Br_train(:,i1,i2,i3,i4,i5);
                    end
                end
            end
        end
    end
    P = ROB_test'* ROB_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5));
    [U,~,Z] = svd(P,'econ');
    Q = U*Z';
    ROB_test = ROB_test*Q;

    %% Step B
    Gamma_K_train = manifold_log(Kr_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5)),Kr_train,'real');
    Gamma_K_test = zeros(N_q,N_q);
    for k = 1:N_q
        for s = 1:N_q
            Gamma_K_test(k,s) = interpn(p(:,:,:,:,:,1),p(:,:,:,:,:,2),p(:,:,:,:,:,3),...
                p(:,:,:,:,:,4),p(:,:,:,:,:,5),...
                reshape(Gamma_K_train(k,s,:,:,:,:,:),size_p1),p1,p2,p3,p4,p5,interpolate_method);
        end
    end
    Kr_test = manifold_exp(Kr_train(:,:,i0(1),i0(2),i0(3),i0(4),i0(5)),Gamma_K_test,'real');

    Gamma_B_train = manifold_log(Br_train(:,i0(1),i0(2),i0(3),i0(4),i0(5)),Br_train,'real');
    Gamma_bias_test = zeros(N_q,1);
    for k = 1:N_q
        Gamma_bias_test(k) = interpn(p(:,:,:,:,:,1),p(:,:,:,:,:,2),p(:,:,:,:,:,3),...
                            p(:,:,:,:,:,4),p(:,:,:,:,:,5),...
                reshape(Gamma_B_train(k,:,:,:,:,:),size_p1),p1,p2,p3,p4,p5,interpolate_method);
    end
    Br_test = manifold_exp(Br_train(:,i0(1),i0(2),i0(3),i0(4),i0(5)), Gamma_bias_test, 'real');
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

function [ROB, Kr,Br] = DMD(X1,X2,r)
    X1_tilde = [X1;ones(1,size(X1,2))];
    [ROB,~,~] = svd(X1,'econ');
    ROB = ROB(:,1:r);
    [V,Sigma,Z] = svd(X1_tilde,'econ');
    V = V(:,1:r+1); Sigma = Sigma(1:r+1,1:r+1); Z = Z(:,1:r+1);
    V1 = V(1:size(X1,1),:);V2 = V(end,:);

    Kr = ROB'*X2*Z/Sigma*V1'*ROB;
    Br = ROB'*X2*Z/Sigma*V2';
    
%     A_tilde = X2*Z/Sigma*V';
%     K = A_tilde(:,1:end-1);
%     B = A_tilde(:,end);
%     Kr = ROB'*K*ROB;
%     Br = ROB'*B;
end

function visualize(T)
    %% domain set up
    x = [167,167,349,349,462,462,260,167,167,258,561,652,652,455,455,341,341,562,652,652,547,269];
    y = [146,295,295,226,226,327,327,417,691,782,782,690,542,542,603,603,508,508,418,146,42,42];
    y = y-10;

    thermalmodel = createpde('thermal','transient');
    R1 =  [3;4;0;800;800;0;0;0;800;800];
    n = length(x);
    P1 = [2;n;x';y'];
    R1 = [R1;zeros(length(P1)-length(R1),1)];
    gd = [R1,P1]; 
    sf = 'R1-P1';
    ns = char('R1','P1');
    ns = ns';
    g = decsg(gd,sf,ns);
    geometryFromEdges(thermalmodel,g);
    thermalProperties(thermalmodel,'ThermalConductivity',1,...
                                   'MassDensity',1,...
                                   'SpecificHeat',1);
    internalHeatSource(thermalmodel,0);

    %% BC
    thermalBC(thermalmodel,'Edge',1,'Temperature',1);
    thermalBC(thermalmodel,'Edge',23,'Temperature',2);
    thermalBC(thermalmodel,'Edge',3:22,'Temperature',3);
    thermalBC(thermalmodel,'Edge',25:26,'Temperature',3);

    %% IC
    thermalIC(thermalmodel,0);

    msh = generateMesh(thermalmodel,'Hmax',50);

    pos1 = [0.08 0.2 0.19 0.7];
    pos2 = [0.41 0.2 0.19 0.7];
    pos3 = [0.74 0.2 0.19 0.7];
    c1_pos = [0.285 0.2 0.01 0.7];
    c2_pos = [0.615 0.2 0.01 0.7];
    c3_pos = [0.945 0.2 0.01 0.7];

    figure
    width = 6.6;     % Width in inches
    height = 2;    % Height in inches
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);

    subplot('Position',pos1);
    pdeplot(thermalmodel,'XYData',T(:,50),'ColorMap','jet')
    xlim([0,800]);
    ylim([0,800]);
    xlabel('$x$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    ylabel('$y$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    hcb = colorbar;
    title(hcb,{'$s$'},'FontUnits','points','interpreter','latex',...
        'FontSize',8);
    set(hcb,'Position',c1_pos,'FontSize',8)

    subplot('Position',pos2)
    pdeplot(thermalmodel,'XYData',T(:,250),'ColorMap','jet')
    xlim([0,800]);
    ylim([0,800]);
    xlabel('$x$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    ylabel('$y$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    hcb = colorbar;
    title(hcb,{'$s$'},'FontUnits','points','interpreter','latex',...
        'FontSize',8);
    set(hcb,'Position',c2_pos,'FontSize',8)

    subplot('Position',pos3)
    pdeplot(thermalmodel,'XYData',T(:,end),'ColorMap','jet')
    xlim([0,800]);
    ylim([0,800]);

    xlabel('$x$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    ylabel('$y$','FontUnits','points','interpreter','latex',...
        'FontSize',10);
    hcb = colorbar;
    title(hcb,{'$s$'},'FontUnits','points','interpreter','latex',...
        'FontSize',8);
    set(hcb,'Position',c3_pos,'FontSize',8);
    

end

