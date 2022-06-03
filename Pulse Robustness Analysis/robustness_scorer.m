function r_score = robustness_scorer(Params)
curdir = 'C:\Users\rmish\OneDrive - Johns Hopkins\Schulman Lab\Topological Sampling Project\5.26.20_Robustness_Analysis';
addpath(curdir,'Model_functions');

Params = [25 25 0 5 0 0 25 0 0 1 0 0 1 1 0 0 -1 0 0 0];
ortho_nodes = 3;

load('k_LHS_mat.mat', 'k_LHS_mat');
k = k_LHS_mat;
%% Model Initialization
time = (0:0.01:3.23)*3600; % seconds
G_top = zeros(length(k), length(time));

for y = 1:length(k) %sample through all k values
    if rem(y,10) == 0
        disp(['Run ' num2str(y)]);
    end
    kpr = k(y,1)*ones(ortho_nodes,1);
    kpc = k(y,2)*ones(ortho_nodes,1);
    kdr = k(y,3)*ones(ortho_nodes,1);
    kdc = k(y,4)*ones(ortho_nodes,1);
    kga = k(y,5)*ones(ortho_nodes,1);
    kgar = k(y,6)*ones(ortho_nodes,1);
    kar = k(y,7)*ones(ortho_nodes,1);
    kgb = k(y,8)*ones(ortho_nodes,1);
    kgbc = k(y,9)*ones(ortho_nodes,1);
    kbc = k(y,10)*ones(ortho_nodes,1);
    kgab = k(y,11)*ones(ortho_nodes,1);
    j = 1;
    rtop = Params(j,end-8:end); %replace param mat with the right value for :>i
    rtop = reshape(rtop,[3 3]);
    [act_vec, blk_vec, CAprod_mat, Rprod_mat, ortho_nodes, ind_genes, init] = gen_mat_converter(rtop);
    rep_vec = act_vec; % vector holding the shared repressors
    ca_vec = blk_vec; % vector holding the shared coactivators
    
    [act_mat] = con_vec2mat(act_vec,ortho_nodes,ind_genes);
    [blk_mat] = con_vec2mat(blk_vec,ortho_nodes,ind_genes);
    [rep_mat] = con_vec2mat(rep_vec,ortho_nodes,ind_genes);
    [ca_mat] = con_vec2mat(ca_vec,ortho_nodes,ind_genes);
    n = ortho_nodes;
    
    
    %% Definining the Parameters from Latin Hypercube Sampling
    %init: figures out which genelets are used (out of 9 possible
    %genelets) 1: used genelet connection 0: unused genelet connection
    G_tot = Params(j, 1:3*n)';
    G_tot(~init) = [];
    % Genelets initial states:  1-ON, 0-BLK
    %if the genelet is ON, remove all blockers. If the genelet is blocked,
    %cause all genelet to be blocked.
    int_vec = Params(j, 3*n+1:4*n); %getting the 1 or 0 from the Parameter list
    init2 = init;
    init2(init == 0)=NaN;
    
    G_int_vec = [init2(1:3).*int_vec(1) init2(4:6).* int_vec(2) init2(7:9).*int_vec(3)];
    G_int_vec(isnan(G_int_vec))=[];
    G_int_vec(G_int_vec == 0) = -1; %no more off genelets - only ON or BLOCKED
    
    dA_tot = 5*(act_mat*G_tot);
    dB_tot = 5*(blk_mat*G_tot); %multiply dB_tot by initial genelet state (binary)
    
    [GdAin,GdBin] = int_genelet_states(G_int_vec,G_tot);
    dBin = zeros(1,ortho_nodes);
    
    for i = 1:3
        z = find(blk_vec == i);
        dBminus = 0;
        for k2 = 1:length(z)+1
            dBminus = dBminus + GdBin(k2);
        end
        if any(dBminus)
            dBin(i) = dB_tot(i) - dBminus; %remove all GdBin from the active blockers
        end
    end
    
    dBin(int_vec == 1) = 0; %take away all blockers from active nodes
    rRin = zeros(1,ortho_nodes);
    dAin = (dA_tot');
    rCAin = zeros(1,ortho_nodes);
    
    int = [rRin,dAin,rCAin,dBin,GdAin,GdBin];
    
    %% Simulating the ODEs
    rates = @(t,x) general_genelet_eqs_topol(t,x,ortho_nodes,ind_genes,dA_tot,dB_tot,G_tot,...
        kpr,kpc,kdr,kdc,kga,kgar,kar,kgb,kgab,kgbc,kbc,...
        act_mat,rep_mat,blk_mat,ca_mat,CAprod_mat,Rprod_mat);
    
    opts = odeset('NonNegative', 1:length(4*ortho_nodes+2*ind_genes)); %#ok<SZARLOG>
    
    [t,x] = ode15s(rates,time,int,opts);
    
    %node 3 start:
    n3s = sum(init(1:6))+1;
    n3e = sum(init);
    %         G3_norm = x(:,4*n+4)/G_tot(4);
    
    G3v = G_tot(n3s:n3e)';
    G3vd = repmat(G3v, length(x), 1);
    
    G3_1 = x(:,4*n+n3s:4*n+n3e);
    
    
    %If the solver breaks, make the terms zeros:
    if length(G3_1) < length(time)
        G3_norm = zeros(length(time),1);
    else
        G3_norm = mean(G3_1./G3vd,2);
    end
    G_top(y,:) = G3_norm;
end

[pcmin, pcmax] = pulsecounter_matrix(G_top,0.1);
r_score = pcmin;
% G_mean = mean(G_top, 1);
% G_var = var(G_top,1);
% avg_var = mean(G_var);
% V = avg_var;
plot(time./3600, G_top');
set(gca, 'FontSize', 35);
% figure;
% shadedErrorBar(time./3600,G_mean,G_var)
% set(gca, 'FontSize', 35);
end